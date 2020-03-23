from setup import *
import cv2
import os
from scipy.signal import medfilt
import traceback

def make_file(tdc, fname):
    run_list, event_list = [], []
    X, Y, Z = [], [], []

    coord0, coord1, coord2, coord3 = [], [], [], []

    if tdc == 24:
        ### Apply cuts bulk24new and singles 24
        for i,B in enumerate(tree):
            # Multiplicity cuts
            passed_mult_cut = False
            if B.nbub == 1 and B.ibub == 1:
                if B.dytranCZT > 0.8 and B.dytranCZT < 1.2:
                    if B.piezo_t0[0] < -0.005 and B.piezo_t0[0] > -0.040:
                        if B.piezo_prehit_filtered_rms[0] < 0.013 and B.piezo_prehit_filtered_rms[2] < 0.0105:
                            if B.piezo_E[39] < 31500 and B.piezo_E[41] < 23500:
                                passed_mult_cut = True

            # Fiducial cuts
            base24 = False
            if B.pset == 30 and abs(B.pts[0]-B.pset) < 1 and abs(B.ts[2] - 16.05) < 0.5 and B.te > 25:
                if B.trigger_main == 0 and B.timestamp > 0 and B.pset > 0 and B.pts[0] > 0 and B.ts[2] > 0:
                    base24 = True

            TazoCuts = False
            Tazo2 = B.TazoProcessed == 2 and B.cDwall < 10 and B.z_theta < -900 and B.cZ > 0
            if B.TazoProcessed != 0 and B.TazoProcessed != 3 and not Tazo2:
                TazoCuts = True

            cutFid24a = False
            if (B.cZ < 0 and B.cDwall > 5) or (B.cZ > 0 and B.cZ < 519 and B.cDwall > 7):
                cutFid24a = True
            if B.cZ > 519 and B.cZ < 529 and B.gFrameDiff < 10 and B.cDwall > 6:
                cutFid24a = True
            if B.cZ > 0 and B.cZ < 519 and B.cDwall < 7 and B.cDwall > 3 and B.z_theta < 0.11 and B.gFrameDiff < 10:
                cutFid24a = True

            passed_fid_cut = False
            if cutFid24a and TazoCuts and base24:
                passed_fid_cut = True

            if passed_mult_cut and passed_fid_cut:
                run_list.append(B.run[:10])
                event_list.append(B.ev)
                X.append(B.X)
                Y.append(B.Y)
                Z.append(B.Z)
                coord0.append([B.vert0, B.hori0])
                coord1.append([B.vert1, B.hori1])
                coord2.append([B.vert2, B.hori2])
                coord3.append([B.vert3, B.hori3])

    else:
        ### Apply cuts bulk33new and singles33
        for i,B in enumerate(tree):
            # Multiplicity cuts
            passed_mult_cut = False
            if B.nbub == 1 and B.ibub == 1:
                if B.dytranCZT > 0.8 and B.dytranCZT < 1.2:
                    if B.piezo_t0[0] < -0.005 and B.piezo_t0[0] > -0.045:
                        if B.piezo_prehit_filtered_rms[0] < 0.013 and B.piezo_prehit_filtered_rms[2] < 0.011:
                            if B.piezo_E[39] < 31500 and B.piezo_E[41] < 27500:
                                passed_mult_cut = True

            # Fiducial cuts
            base33 = False
            if B.pset == 30 and abs(B.pts[0]-B.pset) < 1 and abs(B.ts[2] - 14.05) < 0.5 and B.te > 25:
                if B.trigger_main == 0 and B.timestamp > 0 and B.pset > 0 and B.pts[0] > 0 and B.ts[2] > 0:
                    base33 = True

            TazoCuts = False
            Tazo2 = B.TazoProcessed == 2 and B.cDwall < 10 and B.z_theta < -900 and B.cZ > 0
            if B.TazoProcessed != 0 and B.TazoProcessed != 3 and not Tazo2:
                TazoCuts = True

            cutFid33a = False
            if (B.cZ < 0 and B.cDwall > 10) or (B.cZ >= 0 and B.cZ < 523 and B.cDwall > 11):
                cutFid33a = True
            if B.cZ < 0 and B.cDwall < 10 and B.cDwall > 5 and B.z_theta < 0.12 and B.z_theta > 0:
                cutFid33a = True
            if B.cZ >= 0 and B.cZ < 523 and B.cDwall < 11 and B.cDwall > 4:
                if B.z_theta < 0.12 and B.z_theta_chisq < 20 and B.gFrameDiff < 10:
                    cutFid33a = True

            passed_fid_cut = False
            if cutFid33a and TazoCuts and base33:
                passed_fid_cut = True

            if passed_mult_cut and passed_fid_cut:
                run_list.append(B.run[:10])
                event_list.append(B.ev)
                X.append(B.X)
                Y.append(B.Y)
                Z.append(B.Z)
                coord0.append([B.vert0, B.hori0])
                coord1.append([B.vert1, B.hori1])
                coord2.append([B.vert2, B.hori2])
                coord3.append([B.vert3, B.hori3])


    run_ev_list = list(zip(run_list, event_list))
    run_list = np.array(run_list)
    event_list = np.array(event_list)

    with open(fname, 'w') as f:
        for i, run_ev in enumerate(run_ev_list):
            f.write('%s %d %f %f %f %f %f %f %f %f %f %f %f\n' % (run_ev[0], run_ev[1], coord0[i][0], coord0[i][1],
                    coord1[i][0], coord1[i][1], coord2[i][0], coord2[i][1], coord3[i][0], coord3[i][1], X[i], Y[i], Z[i]))
