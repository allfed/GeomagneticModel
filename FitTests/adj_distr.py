import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bootstrap


adjustments = [
    1.1967969533281275,
    0.8596958507086939,
    0.9364741634508743,
    1.4301531473093996,
    0.8879367120695318,
    1.2306580529581668,
    1.12631848503895,
]

print(np.mean(adjustments), np.std(adjustments), np.ptp(adjustments) / 2)
res = bootstrap((adjustments,), np.mean, confidence_level=0.95, vectorized=False)
l = res.confidence_interval[0]
h = res.confidence_interval[1]
m = (l + h) / 2
d = (h - l) / 2
print(m, d, l, h)
