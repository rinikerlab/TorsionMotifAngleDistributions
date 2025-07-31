#! /bin/env python

from dataRetrieval import GetSmarts, GetHierarchyclass, \
    RetrieveRawDataCsdModified, RetrieveRawDataDash, RetrieveRawDataGNN, RetrieveEnsembleDataVac
from wasserstein import ComputeCramerVonMises, ComputeW1, ComputeW2