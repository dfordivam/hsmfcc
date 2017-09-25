module MFCC where

-- refer http://www.practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/

import Protolude
import Data.Vector.Unboxed
import Numeric.FFT

inputAudioSignalSamplingFreq = 16000 -- 16kHz
frameLengthSeconds = 0.025 -- 25 ms

frameLength = floor (inputAudioSignalSamplingFreq * frameLengthSeconds) -- 400

frameStepSeconds = 0.010 -- 10 ms
frameStepLength = floor (inputAudioSignalSamplingFreq * frameStepLength) -- 160

-- 512 for frameLength of 400
fftSize = val
  where
    powers = map (\x -> 2 ^ x) [1..]
    (val:_) = filter (> frameLength) powers

melFilterBankCount = 40
samplePeriod = 625 -- x100ns for 1 frame

type Frame = Vector Float

data MfccConfig = MfccConfig
  {
  }

computeMfcc :: MfccConfig -> AudioData -> MfccData
computeMfcc =
  let
    frames = createAudioDataFrames
    coefs = extractCoeficients frames

createAudioDataFrames :: AudioData -> Frames

-- number of coefficients per frame 12
-- compute MFCC from windowed data
extractCoeficientsForOneFrame :: Frame -> MfccVector
extractCoeficientsForOneFrame =
  let
    -- Process frame
    p1 = applyHammingWindow . preEmphasize . zMeanFrame

    -- keep only first 257 coefficients
    -- Periodogram based power spectral estimate
    makeFilterBank -- do FFT

    -- Apply DCT to filterbank to make MFCC.
    computeCepstralCoefficients

  -- /* weight cepstrum */
  -- WeightCepstrum(mfcc, para, w);

zMeanFrame :: Frame -> Frame
zMeanFrame frame = fmap (\x -> x - mean) frame
  where mean = (sum frame) / (length frame)

preEmphasize :: Frame -> Float -> Frame
preEmphasize frame preEmph = generate (length frame) f
  where f 0 = (1.0 - preEmph) * (frame ! 0)
        f n = (frame ! n) - ((frame ! (n - 1)) * preEmph)

applyHammingWindow :: Frame -> Frame
applyHammingWindow frame = zipWith (*) frame hammingWindowTable
-- Take the value of banks and then take the log of its power value
-- default
--   number of Filter Banks = 24
--   start freq = 300Hz
--   end freq

makeFilterBank :: Frame -> MfccVector
makeFilterBank frame =

  where
    fOut = fft paddedFrame
    paddedFrame = frameComplex ++ (replicate (fftSize - (length frame)) (0 :+ 0))
    frameComplex = map (\x -> x :+ 0) (V.toList frame)

    pwr = map (\(r :+ i) -> sqrt (r*r + i*i)) fOut

    computeFilterBank
    freq -> bin
    lowerChannelWeights

-- Get the value of power in a bin using triangular filter
filterBin :: Int -> PowerSpecEstimates -> Float
filterBin bin powerVec = sum (zipWith (*) spliceVec weightsVec)

-- (Bin number, bin start freq, bin end freq, weights vector)
type FilterBin = (Int, Int, Int, Vector Float)
type MelBankFilters = Vector FilterBin

melBankFilters :: MelBankFilters
melBankFilters = generate melFilterBankCount f
  where
    f n = (k,loF,hiF,weightsVec)
      where
        k = n + 1 -- 1 to melFilterBankCount
        -- It is the center of previous mel bank
        loF = floor $ mel $ k - 1
        hiF = floor $ mel $ k + 1
        centerF = floor $ mel k
        weightsVec = generate (hiF - loF) g
          where
            -- Triangular
            g m = if freq < centerF
                     then (freq - loF) / (centerF - loF)
                     else (hiF - freq) / (hiF - centerF)
              where
                freq = m + loF

computeFilterBank :: PowerSpecEstimates -> [MelFilterBank]

-- default CepstalCoeficient count = 12
computeCepstralCoefficients :: [MelFilterBank] -> [CepstalCoeficient]
computeCepstralCoefficients melBanks =
  let
    cs = discreetCosineTransform melBanks
  in take 12 cs

-- (gdb) p mhi
-- $35 = 2840.03784

mel :: Int -> Float
mel k = 1127 * log (1 + (k - 1) * fres)
  where fres = inputAudioSignalSamplingFreq / (fftSize * 700)

-- Limit the FFT freq to lower half
nv2 = floor (fftSize/2)

-- Start of first bin freq
mlo = 0
-- End of last bin freq
mhi = mel (nv2 + 1)

-- (gdb) p *mfcc->para
-- {basetype = 7, smp_period = 625, smp_freq = 16000, framesize = 400, frameshift = 160, preEmph = 0.970000029, lifter = 22, fbank_num = 40, delWin = 2, accWin = 2, silFloor = 50, escale = 0.100000001, hipass = -1, lopass = -1,
--   enormal = 0, raw_e = 1, zmeanframe = 1, usepower = 0, vtln_alpha = 1, vtln_upper = -1, vtln_lower = -1, delta = 1, acc = 1, energy = 0, c0 = 0, absesup = 0, cmn = 1, cvn = 1, mfcc_dim = 40, baselen = 40, vecbuflen = 120,
--   veclen = 120, loaded = 0}

-- float Mel(int k, float fres)
-- {
--   return(1127 * log(1 + (k-1) * fres));
-- }

-- Generate table for hamming window.
hammingWindowTable :: Vector Float
hammingWindowTable = generate defaultFrameSize f
  where f n = 0.54 - (0.46 * cos (a * n))
        a = (2.0 * pi)/(defaultFrameSize - 1)
  -- a = 2.0 * PI / (framesize - 1);
  -- for(i=1;i<=framesize;i++) {
  --   /*costbl_hamming[i-1] = 0.54 - 0.46 * cos(2 * PI * (i - 1) / (float)(framesize - 1));*/
  --   w->costbl_hamming[i-1] = 0.54 - 0.46 * cos(a * (i - 1));
  -- }

 -- * Build tables for FFT.
fftTable
  -- for (m = 1; m <= n; m++) {
  --   me = 1 << m;
  --   me1 = me / 2;
  --   w->costbl_fft[m-1] =  cos(PI / me1);
  --   w->sintbl_fft[m-1] = -sin(PI / me1);
  -- }

 -- * Generate table for DCT operation to make mfcc from fbank.
 --  B = PI / fbank_num;
 --  k = 0;
 --  for(i=1;i<=mfcc_dim;i++) {
 --    C = i * B;
 --    for(j=1;j<=fbank_num;j++) {
 --      w->costbl_makemfcc[k] = cos(C * (j - 0.5));
 --      k++;
 --    }
 --  }

 -- * Generate table for weighing cepstrum.
 --  if (lifter > 0) {
 --    a = PI / lifter;
 --    b = lifter / 2.0;
 --    for(i=0;i<mfcc_dim;i++) {
 --      w->sintbl_wcep[i] = 1.0 + b * sin((i+1) * a);
 --    }
 --  } else {
 --    for(i=0;i<mfcc_dim;i++) {
 --      w->sintbl_wcep[i] = 1.0;
 --    }
 --  }

-- void PreEmphasise (float *wave, int framesize, float preEmph)
-- {
--   int i;

--   for(i = framesize; i >= 2; i--)
--     wave[i] -= wave[i - 1] * preEmph;
--   wave[1] *= 1.0 - preEmph;
-- }
