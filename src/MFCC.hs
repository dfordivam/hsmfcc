module MFCC where

-- refer http://www.practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/

import Protolude
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM
import Data.Vector (Vector)
import Numeric.FFT
import Codec.Audio.Wave
import qualified Data.ByteString as BS
import qualified Data.ByteString.Builder as BS
import Control.Monad.ST

import GHC.IO.Handle (hSetPosn, HandlePosn(..))

-- Discreet freq
type Freq = Int

inputAudioSignalSamplingFreq :: Int
inputAudioSignalSamplingFreq = 16000 -- 16kHz

frameLengthSeconds :: Double
frameLengthSeconds = 0.025 -- 25 ms

frameLength :: Int
frameLength = floor ((fromIntegral inputAudioSignalSamplingFreq) * frameLengthSeconds) -- 400

frameStepSeconds :: Double
frameStepSeconds = 0.010 -- 10 ms

frameStepLength :: Int
frameStepLength = floor ((fromIntegral inputAudioSignalSamplingFreq) * frameStepSeconds) -- 160

-- 512 for frameLength of 400
fftSize = val
  where
    powers = map (\x -> 2 ^ x) [1..]
    (val:_) = filter (> frameLength) powers

melFilterBankCount = 40
-- samplePeriod = 625 -- x100ns for 1 frame

type Frame = Vector Double

data MfccConfig = MfccConfig
  {
  }

preEmph = 0.97
processFrame :: Frame -> Frame
processFrame =
  applyHammingWindow . (preEmphasize preEmph) . zMeanFrame

readWaveFileData :: FilePath -> IO ByteString
readWaveFileData fp = do
  -- TODO catch exception
  wd <- readWaveFile fp
  h <- openFile fp ReadMode
  hSetPosn (HandlePosn h (fromIntegral (waveDataOffset wd)))
  bs <- BS.hGetContents h
  BS.writeFile "bs.out" bs
  return bs

getFloat :: Word8 -> Word8 -> Double
getFloat w2 w1 = fromIntegral res2
  where w1shifted :: Word16
        w1shifted = (fromIntegral w1) `shiftL` 8
        res2 :: Int16
        res2 = (fromIntegral res)
        res :: Word16
        res = w1shifted + (fromIntegral w2)

writeFloatData fp = do
  bs <- readWaveFileData fp
  let mels = topAPI bs
      buil = mconcat $ map (\v -> mconcat $
                           map (BS.floatLE . realToFrac) $ V.toList v) mels
  h <- openFile "melfilter.data" WriteMode
  BS.hPutBuilder h buil
-- writeHTKFile fp = do
--   openFile "out.htk" WriteMode
--   bs <- readWaveFileData fp
--   let mels = topAPI bs
--       samples = length mels
--       sampSize = melFilterBankCount * 4 -- 4 Byte float
--       sampPeriod = 625
--       parmKind = + 7
  
topAPI :: ByteString -> [MelFilterBank]
topAPI bs = map (map log) $
  map (makeFilterBank . processFrame) frames
  where
    frames :: [Frame]
    frames =  map f [0 .. (numOfFrames - 1)]
    f :: Int -> Frame
    f i = V.generate frameLength (g i)

    g i j = if index < BS.length bs
               then getFloat (bs `BS.index` index) (bs `BS.index` (index + 1))
               else 0
      where
        index = (2 * i * frameStepLength) + (2 * j)

    numOfFrames = ceiling $ (fromIntegral dataSize) /
      (fromIntegral frameStepLength)
    -- 16 bit of data
    dataSize = floor ((fromIntegral (BS.length bs)) /2)

makeFrames :: ByteString -> IO ([Frame])
makeFrames bs = do
  writeFile "frames.out" (show frames)
  return frames
  where
    frames :: [Frame]
    frames =  map f [0 .. (numOfFrames - 1)]
    f :: Int -> Frame
    f i = V.generate frameLength (g i)

    g i j = if index < BS.length bs
               then getFloat (bs `BS.index` index) (bs `BS.index` (index + 1))
               else 0
      where
        index = (2 * i * frameStepLength) + (2 * j)

    numOfFrames = ceiling $ (fromIntegral dataSize) /
      (fromIntegral frameStepLength)
    -- 16 bit of data
    dataSize = floor ((fromIntegral (BS.length bs)) /2)

zMeanFrame :: Frame -> Frame
zMeanFrame frame = fmap (\x -> x - mean) frame
  where mean = (sum frame) / (fromIntegral $ length frame)

preEmphasize :: Double -> Frame -> Frame
preEmphasize preEmph frame = V.generate (length frame) f
  where f 0 = (1.0 - preEmph) * (frame V.! 0)
        f n = (frame V.! n) - ((frame V.! (n - 1)) * preEmph)

-- preEmphasize
--   for(i = framesize; i >= 2; i--)
--     wave[i] -= wave[i - 1] * preEmph;
--   wave[1] *= 1.0 - preEmph;

applyHammingWindow :: Frame -> Frame
applyHammingWindow frame = V.zipWith (*) frame hammingWindowTable
-- Take the value of banks and then take the log of its power value
-- default
--   number of Filter Banks = 24
--   start freq = 300Hz
--   end freq

type MelFilterBank = Vector Double
type PowerValues = Vector Double

makeFilterBank :: Frame -> MelFilterBank
makeFilterBank frame =
  getBankPower pwr melBankFilters
  where
    fOut = fft paddedFrame
    paddedFrame = frameComplex ++ (replicate (fftSize - (length frame)) (0 :+ 0))
    frameComplex = map (\x -> x :+ 0) (V.toList frame)

    -- Take the first 257 values
    pwr = V.fromList $ take (nv2 + 1) $
      map (\(r :+ i) -> sqrt (r*r + i*i)) fOut

getBankPower :: PowerValues -> MelFilter -> MelFilterBank
getBankPower pwr filt =
  runST $ do
    v <- VM.replicate melFilterBankCount 0.0
    let genV v _ i (b,wt) = do
          VM.modify v (\x -> x + (pwr V.! i) * wt) (b - 1)
          unless ( b == 1) $
            VM.modify v (\x -> x + (pwr V.! i) * (1 - wt)) (b - 2)
    V.ifoldM (genV v) () filt
    V.freeze v

-- Freq -> (bin number, lower weight)
-- For each freq f; it belongs to the lower region of bin k
-- and higher freq region of bin (k - 1)
-- The weight wt is added to the bin k and (1 - wt) to bin (k - 1)
type MelFilter = Vector (Int, Double)

melBankFilters :: MelFilter
melBankFilters = V.generate (nv2 + 1) genF
  where
    -- Start of first bin freq
    mlo = 0
    -- End of last bin freq
    mhi = mel (nv2 + 1)
    -- Mel freq range for 0 to 256
    ms = mhi - mlo

    bankForMelFreq :: Double -> Int
    bankForMelFreq mf = ceiling (((mf - mlo) * (fromIntegral melFilterBankCount)/ms))

    bankMfreq :: Int -> Double
    bankMfreq b = mlo + ((fromIntegral b) * ms) / (fromIntegral melFilterBankCount)

    genF :: Int -> (Int, Double)
    genF 0 = (1, 0)
    genF f = (b,wt)
      where
        --
        b = bankForMelFreq (mel f)

        nextBankMF = bankMfreq b
        firstBankMF = bankMfreq 1
        prevBankMF = bankMfreq (b - 1)
        wt = if b == 1
          -- Special case
          then 1 - (firstBankMF - (mel f)) / (firstBankMF - mlo)
          else 1 - (nextBankMF - (mel f)) / (nextBankMF - prevBankMF)

-- type CepstalCoeficients = Vector Double
-- -- default CepstalCoeficient count = 12
-- computeCepstralCoefficients :: MelFilterBank -> CepstalCoeficients
-- computeCepstralCoefficients melPwrValues =
--   map computeCepCoeff [0 .. (melFilterBankCount - 1)]
--   where
--     computeCepCoeff i =
--       let v = sum $ zipWith (*) melPwrValues (discreetCosineTable ! i)
--       in v * (sqrt (2.0/melFilterBankCount))

-- int i, j, k;
-- k = 0;
-- /* Take DCT */
-- for(i = 0; i < para->mfcc_dim; i++){
--   mfcc[i] = 0.0;
--   for(j = 1; j <= para->fbank_num; j++)
--     mfcc[i] += w->fbank[j] * w->costbl_makemfcc[k++];
--   mfcc[i] *= w->sqrt2var;
-- }

-- (gdb) p mhi
-- $35 = 2840.03784

-- mel freq for given discreet fft freq
mel :: Int -> Double
mel k = 1127 * log (1 + (fromIntegral k) * fres)
  where fres = (fromIntegral inputAudioSignalSamplingFreq) / (fromIntegral (fftSize * 700))

-- discreet fft freq for given mel freq
imel :: Double -> Int
imel mf = floor (r * 700 * (e - 1))
  where e = exp (mf /1127)
        r = (fromIntegral fftSize) /
          (fromIntegral inputAudioSignalSamplingFreq)

-- imelF :: Int -> Freq
-- imelF k = floor $ (imel k) * (fromIntegral (nv2 + 1)) /
--   (imel (melFilterBankCount + 1))

-- Limit the FFT freq to lower half
nv2 = floor ((fromIntegral fftSize)/2)

-- (gdb) p *mfcc->para
-- {basetype = 7, smp_period = 625, smp_freq = 16000, framesize = 400, frameshift = 160, preEmph = 0.970000029, lifter = 22, fbank_num = 40, delWin = 2, accWin = 2, silFloor = 50, escale = 0.100000001, hipass = -1, lopass = -1,
--   enormal = 0, raw_e = 1, zmeanframe = 1, usepower = 0, vtln_alpha = 1, vtln_upper = -1, vtln_lower = -1, delta = 1, acc = 1, energy = 0, c0 = 0, absesup = 0, cmn = 1, cvn = 1, mfcc_dim = 40, baselen = 40, vecbuflen = 120,
--   veclen = 120, loaded = 0}

-- float Mel(int k, float fres)
-- {
--   return(1127 * log(1 + (k-1) * fres));
-- }

-- Generate table for hamming window.
hammingWindowTable :: Vector Double
hammingWindowTable = V.generate frameLength f
  where
    f :: Int -> Double
    f n = 0.54 - (0.46 * cos (a * (fromIntegral n)))

    a = (2.0 * pi)/(fromIntegral (frameLength - 1))
  -- a = 2.0 * PI / (framesize - 1);
  -- for(i=1;i<=framesize;i++) {
  --   w->costbl_hamming[i-1] = 0.54 - 0.46 * cos(a * (i - 1));
  -- }

-- discreetCosineTable :: Vector (Vector Double)
-- discreetCosineTable = generate melFilterBankCount f
--   where
--     b = pi / melFilterBankCount
--     f i = generate melFilterBankCount g
      -- where
      --   -- The for loops start with 1
      --   c = (i + 1) * b
      --   g j = cos (c * (j + 0.5))

--  Generate table for DCT operation to make mfcc from fbank.
--  B = PI / fbank_num;
--  k = 0;
--  for(i=1;i<=mfcc_dim;i++) {
--    C = i * B;
--    for(j=1;j<=fbank_num;j++) {
--      w->costbl_makemfcc[k] = cos(C * (j - 0.5));
--      k++;
--    }
--  }
