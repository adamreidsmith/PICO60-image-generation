from setup import *
import cv2
import os
from scipy.signal import medfilt
import traceback
import matplotlib.pyplot as plt

import make_data_file

'''
This program applies multiplicity and fiduacial cuts to the data to isolate good singles events.
Images of the events are then processed to remove the background and noise, and the bubbles are
isolated.  32x32 pixel images of the bubbles are saved in a directory.  Bubble images are named
"run-event-cam-n.png" where run is the run number, event is the event number in the run, cam is
the camera number, and n is the bubble number in that camera's image.
'''

# Thermodynamic conditions
tdc = 24

# Directory to write the images to
write_dir = './bub_single_%d/' % tdc
# Root file continaing AutoBub output
abubfile = '/home/arsmith/projects/rrg-kenclark/pico/recon_test/devel/30l-16/output/abub3hs_all.root'
# Directory containing raw images of all events
run_path = '/home/arsmith/projects/rrg-kenclark/pico/30l-16-data/'
# Output size of images
out_dim = (32,32)


# Apply fiducial and multiplicity cuts to isolate good single bubble images
if tdc == 24:
    ### Apply cuts bulk24new and singles 24
    run_list, event_list = [], []
    all_runs, all_events = [], []
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

        all_runs.append(B.run[:10])
        all_events.append(B.ev)

else:
    ### Apply cuts bulk33new and singles33
    run_list, event_list = [], []
    all_runs, all_events = [], []
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

        all_runs.append(B.run[:10])
        all_events.append(B.ev)


run_ev_list = list(zip(run_list, event_list))
run_list = np.array(run_list)
event_list = np.array(event_list)

if tdc not in [24,33]:
    run_ev_list = list(zip(all_runs, all_events))
    run_list = np.array(all_runs)
    event_list = np.array(all_events)


fabub = rt.TFile(abubfile, 'read')
abubtree = fabub.Get('T')
### Create a dictionary of relevant autobub information with (run, event) as keys
abub_dict = {}
for B in abubtree:
    run = B.run[:10]
    if (run, B.ev) in run_ev_list:
        if (run, B.ev) in abub_dict:
            abub_dict[(run, B.ev)].append((B.camera, B.frame0, B.hori, B.vert))
        else:
            abub_dict.update({(run, B.ev) : [(B.camera, B.frame0, B.hori, B.vert)]})
            

image_read_identifier = cv2.IMREAD_GRAYSCALE

def get_run_init_images(run_num):
    #Get first two images in each event of a specific run for each cam
    im_list = [[],[],[],[]]
    for event in os.listdir(run_path + run_num):
        try:
            float(event)
            for cam in range(4):
                im1_path = run_path + run_num + '/' + event + '/Images/cam%d_image30.png' % cam
                im2_path = run_path + run_num + '/' + event + '/Images/cam%d_image31.png' % cam
                im_list[cam].append(cv2.imread(im1_path, image_read_identifier))
                im_list[cam].append(cv2.imread(im2_path, image_read_identifier))
        except ValueError:
            pass
    return im_list


def im_mean_stdev(cv_im_list):
    #Compute the mean and standard deviation of each pixel in a list of cv2 images
    N = len(cv_im_list)
    mean = np.mean(cv_im_list, axis=0)
    stdev = np.std(cv_im_list, axis=0)
    return mean, stdev


run_means, run_stdevs = {}, {}
def run_means_stdevs(run_num):
    #Return the run image means and stevs of they have already been computed
    if run_num in run_means:
        return run_means[run_num], run_stdevs[run_num]
    
    #Otherwise, compute them
    init_images = get_run_init_images(run_num)

    means, stdevs = [], []
    for cam in range(4):
        mean, stdev = im_mean_stdev(init_images[cam])
        means.append(mean)
        stdevs.append(stdev)

    #Update the dictionaries
    run_means.update( {run_num : means} )
    run_stdevs.update( {run_num : stdevs} )

    return means, stdevs


def get_gen_frames(run_num, event):
    #Get the first frame in which a bubble is seen
    event_data = abub_dict[(run_num, event)]
    gen_frame_nums = [bub_data[1] for bub_data in event_data]
    max_frame_num = max(gen_frame_nums)  #Latest genesis frame from all cameras
    gen_frames = []
    for cam in range(4):
        frame = gen_frame_nums[cam] if gen_frame_nums[cam] > 29 else max_frame_num
        if cam == 2:
            frame += 1
        im_path = run_path + run_num + '/' + str(event) + '/' + '/Images/cam%d_image%d.png' % (cam, frame)
        gen_frames.append(cv2.imread(im_path, image_read_identifier))
    return gen_frames


def isolate_bub(trig_frame, mean, stdev):
    #Isolate bubble in the trigger/genesis frame as per Pitam's algorithm

    #Pixels which are at least 6 sigma above the noise
    F6sigma = np.abs(trig_frame - mean) - 6*stdev

    #Apply a median blur to the image
    im_blur = medfilt(F6sigma, kernel_size=3)

    #Set every pixel with intensity <= 3 to zero
    im_thres = np.where(im_blur <= 3, 0, im_blur)

    #Apply thresholding again with Otsu's binarization
    thres, im_iso = cv2.threshold(im_thres.astype('uint8'), thresh=2, maxval=255, type=cv2.THRESH_OTSU)

    return im_iso


def get_bub_locs(run_num, event, cam):
    all_bubs = abub_dict[(run_num, event)]
    bub_locs = []

    for bub in all_bubs:
        if bub[0] == cam:
            bub_locs.append((round(bub[2]), round(bub[3])))
    return bub_locs


'''
#Remove background and noise to analyze images and save them in a directory
os.makedirs('./full_single_%d/' % tdc, exist_ok=True)
im_done = os.listdir('full_single_%d/' % tdc)
for run, event in run_ev_list:

    done = False
    for im in im_done:
        if run + '_' + str(event) in im:
            done = True
    if done:
        print(run + '_' + str(event) + ' done')
        continue
    
    try:
        print('Analyzing run %s event %d' % (run, event))
    
        #Get the pixel mean and standard deviation of the run
        means, stdevs = run_means_stdevs(run)
        genesis_frames = get_gen_frames(run, event)

        for cam in range(4):
            #Isolate the bubble in the event
            bub_iso = isolate_bub(genesis_frames[cam], means[cam], stdevs[cam])

            #Continue if image contains no bubble (i.e. all pixels are black)
            if np.sum(bub_iso) <= 1:
                continue

            im_name = 'full_single_%d/%s-%d-%d.png' % (tdc, run, event, cam)
            cv2.imwrite(im_name, bub_iso)

    except Exception:
        traceback.print_exc()
        print('Analysis of run %s event %d failed' % (run, event))
'''

#Isolate bubbles in all singles events and save 32x32 pixel images of the bubbles
failed = []
nfailed = 0
os.makedirs(write_dir, exist_ok=True)
im_done = os.listdir(write_dir)
for run, event in run_ev_list:

    done = False
    for im in im_done:
        if run + '-' + str(event) in im:
            done = True
    if done:
        print(run + '-' + str(event) + ' done')
        continue

    try:
        print('Analyzing run %s event %d' % (run, event))
    
        #Get the pixel mean and standard deviation of the run
        means, stdevs = run_means_stdevs(run)
        genesis_frames = get_gen_frames(run, event)

        for cam in range(4):
            #Isolate the bubble in the event
            bub_iso = isolate_bub(genesis_frames[cam], means[cam], stdevs[cam])
            
            #Get locations of all bubbles in the event
            box_locs = get_bub_locs(run, event, cam)
            
            n = 0
            for box_loc in box_locs:
                #Continue if bubble position not known
                if box_loc == (0,0):
                    continue

                #Crop the image around the bubble
                miny, maxy = box_loc[1] - out_dim[1]//2, box_loc[1] + out_dim[1]//2
                minx, maxx = box_loc[0] - out_dim[0]//2, box_loc[0] + out_dim[0]//2
                out_im = bub_iso[miny:maxy, minx:maxx]

                #Continue if image contains no bubble (i.e. all pixels are black)
                if np.sum(out_im) <= 1:
                    continue

                #Write the image to the specified directory
                bub_name = write_dir + '%s-%d-%d-%d.png' % (run, event, cam, n)
                print('Writing bubble to ' + bub_name)
                cv2.imwrite(bub_name, out_im)
                n += 1

    except Exception:
        print('Analysis of run %s event %d failed' % (run, event))
        traceback.print_exc()
        if (run, event) not in failed:
            failed.append((run, event))
        nfailed += 1

print('Anlysis complete. Events failed: ', failed)
print('Number of failed events: ', nfailed)

print('Making data file')
make_data_file.make_file(tdc, 'bubble_data_%d.txt' % tdc)
    
