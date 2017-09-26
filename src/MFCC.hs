module MFCC where

-- refer http://www.practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/

import Protolude
import qualified Data.Vector as V
import Data.Vector (Vector)
import Numeric.FFT
import Codec.Audio.Wave
import qualified Data.ByteString as BS

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
getFloat w1 w2 = fromIntegral (w1shifted + (fromIntegral w2))
  where w1shifted :: Word16
        w1shifted = (fromIntegral w1) `shiftL` 8

topAPI :: ByteString -> [MelFilterBank]
topAPI bs = map (makeFilterBank . processFrame) frames
  where
    frames :: [Frame]
    frames =  map f [0 .. (numOfFrames - 1)]
    f :: Int -> Frame
    f i = V.generate frameLength $ \j ->
      getFloat (bs `BS.index` ((2 * i * frameStepLength) + (2 * j)))
        (bs `BS.index` ((2 * i * frameStepLength) + (2 * j) + 1))


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
  map (getPowerInBank pwr) melBankFilters
  where
    fOut = fft paddedFrame
    paddedFrame = frameComplex ++ (replicate (fftSize - (length frame)) (0 :+ 0))
    frameComplex = map (\x -> x :+ 0) (V.toList frame)

    -- Take the first 257 values
    pwr = V.fromList $ take (nv2 + 1) $
      map (\(r :+ i) -> sqrt (r*r + i*i)) fOut

getPowerInBank :: PowerValues -> FilterBin -> Double
getPowerInBank pwr (k,loF,hiF,weightsVec) =
  sum $ V.zipWith (*) weightsVec pwrBin
  where pwrBin = V.slice loF hiF pwr

-- (Bin number, bin start freq, bin end freq, weights vector)
type FilterBin = (Int, Int, Int, Vector Double)
type MelBankFilters = Vector FilterBin

melBankFilters :: MelBankFilters
melBankFilters = V.generate melFilterBankCount f
  where
    -- Start of first bin freq
    mlo = mel 0
    -- End of last bin freq
    mhi = mel (nv2 + 1)
    -- Step 0 to 256
    ms = mhi - mlo

    -- Mel freq for bank k
    bankMfreq k =
      ((fromIntegral k) * ms)/(fromIntegral melFilterBankCount)

    f n = (k,loF,hiF,weightsVec)
      where
        k = n + 1 -- 1 to melFilterBankCount
        -- Mel freq for this bank

        centerMF = bankMfreq k
        loMF = bankMfreq $ k - 1
        hiMF = bankMfreq $ k + 1
        -- It is the center of previous mel bank
        loF = imel loMF
        hiF = imel hiMF
        centerF = imel centerMF
        weightsVec = V.generate (hiF - loF + 1) g
          where
            -- Triangular
            g m = if freq < centerF
                     then 1 - (c - (mel freq)) /
                          (c - l)
                     else (h - (mel freq)) /
                          (h - c)
              where
                freq = m + loF
                c = mel centerF
                l = mel loF
                h = mel hiF

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
