module Main where


log2 :: Integral a => a -> a
log2 n = ceiling (log (fromIntegral n) / log 2)


log3 :: Integral a => a -> a
log3 n = ceiling (log (fromIntegral n) / log 3)


log4 :: Integral a => a -> a
log4 n = ceiling (log (fromIntegral n) / log 4)


cost_crypto2020 :: Integral a => a -> a
cost_crypto2020 b = 3*b + b*b


cost_new :: Integral a => a -> a
cost_new b = 3*b*log2 b + 2*b


cost_new_3 :: Integral a => a -> a
cost_new_3 b = b*(log3 b + 1) + b * 2 * log3 b + b*log3 b + b


cost_new_4 :: Integral a => a -> a
cost_new_4 b = b*(log4 b + 1) + b * 3 * log4 b + b*log4 b + b


cost_comparison :: Integral a => a -> (a, a)
cost_comparison b = (cost_crypto2020 b, cost_new b)


improvement :: Integral a => a -> Double
improvement b = fromIntegral (cost_crypto2020 b) / fromIntegral (cost_new b)
