# source: https://github.com/Xinglab/rmats-turbo-tutorial/tree/main/scripts

import math
import numpy as np

from class_exon import *
from class_fromGTF import *
from extract_PSI_count import *
from get_novel_ID import *
from rmats_filtering import *


def read_rMATS(fn, readCov, minPSI, maxPSI, sigFDR, sigDeltaPSI, bgFDR, bgWithinGroupDeltaPSI):
    exon, event_type = get_exon_class(fn)
    filtered_event_list = []
    up_event_list = []
    dn_event_list = []
    bg_event_list = []
    with open(fn, "r") as f:
        header = f.readline()
        for line in f:
            x = exon(line)
            if (
                x.averageCountSample1 >= readCov
                and x.averageCountSample2 >= readCov
                and np.nanquantile(x.IncLevel1 + x.IncLevel2, 0.05) <= maxPSI
                and np.nanquantile(x.IncLevel1 + x.IncLevel2, 0.95) >= minPSI
                and len(x.chrom) <= 5
            ):
                filtered_event_list.append(x)
                if x.FDR <= sigFDR:
                    if x.IncLevelDifference >= sigDeltaPSI:
                        dn_event_list.append(x)
                    elif x.IncLevelDifference <= -sigDeltaPSI:
                        up_event_list.append(x)
                elif (
                    x.FDR >= bgFDR
                    and max(x.IncLevel1) - min(x.IncLevel1) <= bgWithinGroupDeltaPSI
                    and max(x.IncLevel2) - min(x.IncLevel2) <= bgWithinGroupDeltaPSI
                ):
                    bg_event_list.append(x)
            # considering cases when --b2 is not provided.
            elif (
                x.averageCountSample1 >= readCov
                and math.isnan(x.averageCountSample2)
                and np.nanquantile(x.IncLevel1 + x.IncLevel2, 0.05) <= maxPSI
                and np.nanquantile(x.IncLevel1 + x.IncLevel2, 0.95) >= minPSI
                and len(x.chrom) <= 5
            ):
                filtered_event_list.append(x)
    event_dict = {
        "upregulated": up_event_list,
        "downregulated": dn_event_list,
        "background": bg_event_list,
        "filtered": filtered_event_list,
    }
    return header, event_dict
