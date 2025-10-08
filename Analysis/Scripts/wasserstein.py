import numpy as np
from scipy.spatial.distance import pdist, squareform, jensenshannon
from scipy.optimize import linprog
import ot
from scipy.stats import cramervonmises_2samp

def ComputeCramerVonMises(dat1, dat2, radian):
    if not radian:
        dat1 = np.deg2rad(dat1)
        dat2 = np.deg2rad(dat2)
    dat1 = np.random.choice(dat1, 1000)
    dat2 = np.random.choice(dat2, 1000)
    symmetrizedDat1 = np.concatenate((dat1, 2*np.pi-dat1))
    symmetrizedDat2 = np.concatenate((dat2, 2*np.pi-dat2))
    globalMin = min(np.min(symmetrizedDat1), np.min(symmetrizedDat2))
    globalMax = max(np.max(symmetrizedDat1), np.max(symmetrizedDat2))
    symmetrizedDat1 = (symmetrizedDat1 - globalMin) / (globalMax - globalMin)
    symmetrizedDat2 = (symmetrizedDat2 - globalMin) / (globalMax - globalMin)
    cm = cramervonmises_2samp(symmetrizedDat1, symmetrizedDat2)
    omega = (cm.statistic/len(symmetrizedDat1))**0.5
    return omega

def ComputeW1(dat1, dat2, radian):
    if not radian:
        dat1 = np.deg2rad(dat1)
        dat2 = np.deg2rad(dat2)
    symmetrizedDat1 = np.concatenate((dat1, 2*np.pi-dat1))
    symmetrizedDat2 = np.concatenate((dat2, 2*np.pi-dat2))
    globalMin = min(np.min(symmetrizedDat1), np.min(symmetrizedDat2))
    globalMax = max(np.max(symmetrizedDat1), np.max(symmetrizedDat2))
    symmetrizedDat1 = (symmetrizedDat1 - globalMin) / (globalMax - globalMin)
    symmetrizedDat2 = (symmetrizedDat2 - globalMin) / (globalMax - globalMin)
    bins = np.linspace(0, 1+0.013, 74)
    hist1, _ = np.histogram(symmetrizedDat1, bins=bins, density=False)
    hist1 = hist1 / hist1.sum()
    hist2, _ = np.histogram(symmetrizedDat2, bins=bins, density=False)
    hist2 = hist2 / hist2.sum()
    bin_centers = (bins[:-1] + bins[1:]) / 2
    C = ot.dist(bin_centers.reshape(-1, 1), bin_centers.reshape(-1, 1), metric='euclidean')
    d_emd2 = ot.emd2(hist1, hist2, C)
    return np.sqrt(d_emd2)

def ComputeW2(dat1, dat2, radian):
    if not radian:
        dat1 = np.deg2rad(dat1)
        dat2 = np.deg2rad(dat2)
    symmetrizedDat1 = np.concatenate((dat1, 2*np.pi-dat1))
    symmetrizedDat2 = np.concatenate((dat2, 2*np.pi-dat2))
    globalMin = min(np.min(symmetrizedDat1), np.min(symmetrizedDat2))
    globalMax = max(np.max(symmetrizedDat1), np.max(symmetrizedDat2))
    symmetrizedDat1 = (symmetrizedDat1 - globalMin) / (globalMax - globalMin)
    symmetrizedDat2 = (symmetrizedDat2 - globalMin) / (globalMax - globalMin)
    bins = np.linspace(0, 1+0.013, 74)
    hist1, _ = np.histogram(symmetrizedDat1, bins=bins, density=False)
    hist1 = hist1 / hist1.sum()
    hist2, _ = np.histogram(symmetrizedDat2, bins=bins, density=False)
    hist2 = hist2 / hist2.sum()
    bin_centers = (bins[:-1] + bins[1:]) / 2
    C = ot.dist(bin_centers.reshape(-1, 1), bin_centers.reshape(-1, 1), metric='sqeuclidean')
    d_emd2 = ot.emd2(hist1, hist2, C)
    return np.sqrt(d_emd2)