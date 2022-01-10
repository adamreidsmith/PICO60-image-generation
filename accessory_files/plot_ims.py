import pickle
import numpy as np
import matplotlib.pyplot as plt
import cv2
#import ks2d

data24_file = './data/pix_pos_rad_dict24.pkl'
data33_file = './data/pix_pos_rad_dict33.pkl'
fake_file = './artificial_events/singles_new_scale_mar21/data_with_rad_dict.pkl'
run_type24_file = './data/run_type24.txt'
run_type33_file = './data/run_type33.txt'

with open(data24_file, 'rb') as f, open(data33_file, 'rb') as g, open(fake_file, 'rb') as h:
    data24 = pickle.load(f)
    data33 = pickle.load(g)
    fake_data = pickle.load(h)
    
fake_pix, fake_pos, fake_rad = [[], [], [], []], [[], [], [], []], [[], [], [], []]
for cam in range(4):
    for key in fake_data:
        if key[-1] == cam and sum(fake_data[key][1][cam]) > 0:
            fake_pix[cam].append(list(fake_data[key][1][cam]))
            fake_pos[cam].append(fake_data[key][0])
            fake_rad[cam].append(fake_data[key][-1])
    fake_pix[cam] = np.array(fake_pix[cam])
    fake_pos[cam] = np.array(fake_pos[cam])
    fake_rad[cam] = np.array(fake_rad[cam])

real_pix24, real_pos24, real_rad24, run_ev24 = [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]
for cam in range(4):
    for key in data24:
        if key[-1] == cam and sum(data24[key][0][cam]) > 0:
            real_pix24[cam].append(list(data24[key][0][cam]))
            real_pos24[cam].append(data24[key][1]/10)  # Position originally in mm
            real_rad24[cam].append(data24[key][-1])
            run_ev24[cam].append(key[:2])
    real_pix24[cam] = np.array(real_pix24[cam])
    real_pos24[cam] = np.array(real_pos24[cam])
    real_rad24[cam] = np.array(real_rad24[cam])

real_pix33, real_pos33, real_rad33, run_ev33 = [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]
for cam in range(4):
    for key in data33:
        if key[-1] == cam and sum(data33[key][0][cam]) > 0:
            real_pix33[cam].append(list(data33[key][0][cam]))
            real_pos33[cam].append(data33[key][1]/10)  # Position originally in mm
            real_rad33[cam].append(data33[key][-1])
            run_ev33[cam].append(key[:2])
    real_pix33[cam] = np.array(real_pix33[cam])
    real_pos33[cam] = np.array(real_pos33[cam])
    real_rad33[cam] = np.array(real_rad33[cam])

# Concatenate all pixel coordinates from real data
real_pix, real_pos, real_rad, run_ev = [], [], [], []
for i in range(4):
    real_pix.append(np.array([*real_pix24[i], *real_pix33[i]]))
    real_pos.append(np.array([*real_pos24[i], *real_pos33[i]]))
    real_rad.append(np.array([*real_rad24[i], *real_rad33[i]]))
    run_ev.append([*run_ev24[i], *run_ev33[i]])


with open(run_type24_file, 'r') as f, open(run_type33_file, 'r') as g:
    lines_24 = f.readlines()
    lines_33 = g.readlines()

run_type24_dict = {}
for line in lines_24:
    line = line.split()
    run_type24_dict.update( { (line[0], int(line[1])) : int(line[2]) } )
    
run_type33_dict = {}
for line in lines_33:
    line = line.split()
    run_type33_dict.update( { (line[0], int(line[1])) : int(line[2]) } )
    
run_type_dict = {**run_type24_dict, **run_type33_dict}

source_pos_in_tube = {0:-1., 10:-1., 100:-1., 14:40., 15:60., 16:80., 21:100., 22:45., 23:120., 24:10., 32:45., 41:45., 99:-1, 2:146.05}
tube_pos = np.array([23.80, -23.80, -79.62])

weights = [[], [], [], []]
for cam in range(4):
    for i, r_e in enumerate(run_ev[cam]):
        try:
            if run_type_dict[r_e] < 0:
                raise KeyError
            source_pos = tube_pos.copy()
            source_pos[2] += source_pos_in_tube[run_type_dict[r_e]]
            vec = real_pos[cam][i] - source_pos
            normalization_factor = np.dot(vec, vec)   # NORMALIZING BY R^2
        except KeyError:
            vec = np.array([0,0,(52.79 - 6.59)/2]) - (tube_pos + np.array([0, 0, 146.05/2]))
            normalization_factor = np.dot(vec, vec)   # NORMALIZING BY R^2
        weights[cam].append(normalization_factor)
    weights[cam] = np.array(weights[cam])
    

def PCC(h1, h2):
    h1 /= np.sum(h1)
    h2 /= np.sum(h2)
    N = len(h1)
    M = len(h1[0])
    num = 0
    den1, den2 = 0, 0
    m1, m2 = h1.mean(), h2.mean()
    for i in range(N):
        for j in range(M):
            num += (h1[i][j] - m1)*(h2[i][j] - m2)
            den1 += (h1[i][j] - m1)**2
            den2 += (h2[i][j] - m2)**2
    return num/(np.sqrt(den1*den2))
            

def plot_hist(pix1, pix2, weights1=None, weights2=None, pix_per_bin1=50, pix_per_bin2=50):
    binsx1 = np.linspace(0,1700,1700//pix_per_bin1)
    binsy1 = np.linspace(0,1088,1088//pix_per_bin1)
    
    binsx2 = np.linspace(0,1700,1700//pix_per_bin2)
    binsy2 = np.linspace(0,1088,1088//pix_per_bin2)
    
    #fig, ((ax11, ax12),(ax21, ax22),(ax31, ax32),(ax41, ax42)) = plt.subplots(4, 2)
    fig, axes = plt.subplots(4, 2, figsize=(14, 14))
    
    cm = plt.cm.viridis

    for cam in range(4):
        (ax1, ax2) = axes[cam]

        if weights1 is not None:
            im1 = ax1.hist2d(pix1[cam][:,0], pix1[cam][:,1], bins=[binsx1, binsy1], cmap=cm)

            weights1[cam] *= np.sum(im1[0])/np.sum(weights1[cam])
            
            z1, xedges, yedges = np.histogram2d(pix1[cam][:,0], pix1[cam][:,1], bins=[binsx1, binsy1], weights=weights1[cam])
            im1 = ax1.pcolormesh(xedges, yedges, z1.transpose(), vmin=0, vmax=52)
            
            #im1 = ax1.hist2d(pix1[cam][:,0], pix1[cam][:,1], weights=weights1[cam], bins=[binsx1, binsy1], cmap=cm)
        else:
            z1, xedges, yedges = np.histogram2d(pix1[cam][:,0], pix1[cam][:,1], bins=[binsx1, binsy1])
            im1 = ax1.pcolormesh(xedges, yedges, z1.transpose(), vmin=0, vmax=52)
        ax1.invert_yaxis()
        
        z2, xedges, yedges = np.histogram2d(pix2[cam][:,0], pix2[cam][:,1], bins=[binsx1, binsy1])
        #im2 = ax2.hist2d(pix2[cam][:,0], pix2[cam][:,1], bins=[binsx1, binsy1], cmap=cm)
        
        r = PCC(z1, z2)
    
        if weights2 is not None:
            z2, xedges, yedges = np.histogram2d(pix2[cam][:,0], pix2[cam][:,1], bins=[binsx2, binsy2], weights=weights2[cam])
            im2 = ax1.pcolormesh(xedges, yedges, z2.transpose(), vmin=0, vmax=160)
            #im2 = ax2.hist2d(pix2[cam][:,0], pix2[cam][:,1], weights=weights2[cam], bins=[binsx2, binsy2], cmap=cm)
        else:
            z2, xedges, yedges = np.histogram2d(pix2[cam][:,0], pix2[cam][:,1], bins=[binsx2, binsy2])
            im2 = ax2.pcolormesh(xedges, yedges, z2.transpose(), vmin=0, vmax=160)
            #im2 = ax2.hist2d(pix2[cam][:,0], pix2[cam][:,1], bins=[binsx2, binsy2], cmap=cm)
        ax2.invert_yaxis()

#        ax1.tick_params(labelbottom=False)
#        ax2.tick_params(labelbottom=False)
#        ax1.tick_params(labelleft=False)
#        ax2.tick_params(labelleft=False)

        cax1 = plt.colorbar(im1, ax=ax1)
        cax1.ax.set_visible(False)
        cax2 = plt.colorbar(im2, ax=ax2)
        cax2.ax.set_visible(False)

        ax2.text(1400,1060,'$r=%.2f$' % r, c='w', size=12)
        
        if cam == 0:
            ax1.set_title('PICO-60 images', size=14)
            ax2.set_title('Artificial images', size=14)
            
    fig.text(0.5, 0.02, 'Horizontal coordinate [pixels]', ha='center', size=14)
    fig.text(0.025, 0.5, 'Vertical coordinate [pixels]', va='center', rotation='vertical', size=14)

    fig.text(0.045, 0.88, 'Camera 0', va='center', rotation='vertical', size=12)
    fig.text(0.045, 0.64, 'Camera 1', va='center', rotation='vertical', size=12)
    fig.text(0.045, 0.4, 'Camera 2', va='center', rotation='vertical', size=12)
    fig.text(0.045, 0.16, 'Camera 3', va='center', rotation='vertical', size=12)
    
#    fig.text(0.5, 0.92, 'Distribution of Bubbles in Artificial and Real PICO-60 Images', ha='center', size=16)
    
    cb_ax1 = fig.add_axes([0.463, 0.058, 0.01, 0.915])
    cbar = fig.colorbar(im1, cax=cb_ax1)
    
    cb_ax2 = fig.add_axes([0.907, 0.058, 0.01, 0.915])
    cbar = fig.colorbar(im2, cax=cb_ax2)

#    cbar.set_ticks([])

    #plt.colorbar()
    plt.tight_layout(rect=[0.05, 0.03, 1, 1], w_pad=3)
    #plt.show()

plot_hist(pix1=real_pix, pix2=fake_pix, weights1=weights, pix_per_bin1=90, pix_per_bin2=40)
plt.savefig('bubble_distribution_cbn.png')

cam_positions = np.array(
                [[-19.0500-0.8351, -32.9956-1.4465, 40.8681],
                [-19.0500-0.7743, -32.9956-1.3411, 13.0119],
                [19.0500+1.2247, -32.9956-2.1212, 40.9293],
                [19.0500+1.6113, -32.9956-2.7908, 12.8003]]
                )
def dist_to_camera(p, cam):
    return np.linalg.norm(cam_positions[cam] - p)

#pts = np.array([[500, 250],
#       [1700//2, 250],
#       [1700-500, 250],
#       [500, 1088//2],
#       [1700//2, 1088//2],
#       [1700-500, 1088//2],
#       [500, 1088-250],
#       [1700//2, 1088-250],
#       [1700-500, 1088-250]]
#       )

#pts = np.array([[600, 400],
#                [1100, 400],
#                [600, 1088-400],
#                [1100, 1088-400]])

#pts = np.array([[650, 1088//2], [1700-650, 1088//2]])
#
#s=[]
#b=[]
#for pt in pts:
#    fig, ax = plt.subplots(4, 2, figsize=(12,16))
#    for cam in range(4):
#        ax1, ax2 = ax[cam]
#        dist_from_cam_real, radius_real = [], []
#        for i, coord in enumerate(real_pix[cam]):
#            if np.linalg.norm(pt - coord) < 200.:
#                dist_from_cam_real.append(dist_to_camera(real_pos[cam][i], cam))
#                radius_real.append(real_rad[cam][i])
#
#        dist_from_cam_fake, radius_fake = [], []
#        for i, coord in enumerate(fake_pix[cam]):
#            if np.linalg.norm(pt - coord) < 200.:
#                dist_from_cam_fake.append(dist_to_camera(fake_pos[cam][i], cam))
#                radius_fake.append(fake_rad[cam][i])
#
#        dist_from_cam_real, dist_from_cam_fake = np.array(dist_from_cam_real), np.array(dist_from_cam_fake)
#        radius_real, radius_fake = np.array(radius_real), np.array(radius_fake)
#
#        ax1.scatter(dist_from_cam_real, radius_real, s=2)
#        ax2.scatter(dist_from_cam_fake, radius_fake, s=2)
#
##        ax1.set_title('Camera %i PICO data' % cam)
##        ax2.set_title('Camera %i artificial data' % cam)
#
#        ax1.set_ylim((0,6))
#        ax2.set_ylim((0,6))
#
##        ax1.set_xlabel('Distance from camera [cm]')
##        ax1.set_ylabel('Bubble radius [pixels]')
##        ax2.set_xlabel('Distance from camera [cm]')
##        ax2.set_ylabel('Bubble radius [pixels]')
#
#        lim = 50
#        real_ind = np.where(radius_real < lim)
#        fake_ind = np.where(radius_fake < lim)
#
#        fit1, cov1 = np.polyfit(dist_from_cam_real[real_ind], radius_real[real_ind], 1, cov=True)
#        fit2, cov2 = np.polyfit(dist_from_cam_fake[fake_ind], radius_fake[fake_ind], 1, cov=True)
#        cov1 = np.sqrt(np.diag(cov1))
#        cov2 = np.sqrt(np.diag(cov2))
#
#        x = np.linspace(25, 60, 200)
#        y1 = fit1[0]*x + fit1[1]
#        y2 = fit2[0]*x + fit2[1]
#
#        ax1.plot(x, y1, c='r', label='$r=(%.3f\pm%.3f)d+(%.1f\pm%.1f)$' % (fit1[0], cov1[0], fit1[1], cov1[1]))
#        ax2.plot(x, y2, c='r', label='$r=(%.3f\pm%.3f)d+(%.2f\pm%.2f)$' % (fit2[0], cov2[0], fit2[1], cov2[1]))
#
#        ax1.legend()
#        ax2.legend()
#
#        s.append(fit1[0])
#        b.append(fit1[1])
#
#        if cam == 0:
#            ax1.set_title('PICO-60 images')
#            ax2.set_title('Artificial images')
#
#    fig.text(0.5, 0.07, 'Distance from camera $d$ [cm]', ha='center', size=12)
#    fig.text(0.045, 0.5, 'Bubble radius $r$ [pixels]', va='center', rotation='vertical', size=12)
#
#    fig.text(0.07, 0.8, 'Camera 0', va='center', rotation='vertical')
#    fig.text(0.07, 0.6, 'Camera 1', va='center', rotation='vertical')
#    fig.text(0.07, 0.4, 'Camera 2', va='center', rotation='vertical')
#    fig.text(0.07, 0.2, 'Camera 3', va='center', rotation='vertical')
#
#    fig.text(0.5, 0.92, 'Radius of bubbles at pixel coordinate [%i, %i]' % (pt[0], pt[1]), ha='center', size=16)
#
#    #plt.show()
#    #plt.tight_layout()
#    plt.savefig('radius_at_(%i,%i).png' % (pt[0], pt[1]))
#
#print(sum(s)/len(s))
#print(sum(b)/len(b))




#s, b = [], []
fig, ax = plt.subplots(4, 2, figsize=(12,16))
w = 150
for cam in range(4):
    ax1, ax2 = ax[cam]
    dist_from_cam_real, radius_real = [], []
    for i, coord in enumerate(real_pix[cam]):
        if 1088//2-w < coord[1] < 1088//2+w:
            if (cam == 1 or cam == 3) and coord[0] > 1250:
                continue
            dist_from_cam_real.append(dist_to_camera(real_pos[cam][i], cam))
            radius_real.append(real_rad[cam][i])

    dist_from_cam_fake, radius_fake = [], []
    for i, coord in enumerate(fake_pix[cam]):
        if 1088//2-w < coord[1] < 1088//2+w:
            if (cam == 1 or cam == 3) and coord[0] > 1250:
                continue
            dist_from_cam_fake.append(dist_to_camera(fake_pos[cam][i], cam))
            radius_fake.append(fake_rad[cam][i])
    
    dist_from_cam_real, dist_from_cam_fake = np.array(dist_from_cam_real), np.array(dist_from_cam_fake)
    radius_real, radius_fake = np.array(radius_real), np.array(radius_fake)

    ax1.scatter(dist_from_cam_real, radius_real, s=2)
    ax2.scatter(dist_from_cam_fake, radius_fake, s=2)
    
#        ax1.set_title('Camera %i PICO data' % cam)
#        ax2.set_title('Camera %i artificial data' % cam)
    
    ax1.set_ylim((0,6))
    ax2.set_ylim((0,6))
    
#        ax1.set_xlabel('Distance from camera [cm]')
#        ax1.set_ylabel('Bubble radius [pixels]')
#        ax2.set_xlabel('Distance from camera [cm]')
#        ax2.set_ylabel('Bubble radius [pixels]')
    
    lim = 6
    real_ind = np.where(radius_real < lim)
    fake_ind = np.where(radius_fake < lim)
    
    fit1, cov1 = np.polyfit(dist_from_cam_real[real_ind], radius_real[real_ind], 1, cov=True)
    fit2, cov2 = np.polyfit(dist_from_cam_fake[fake_ind], radius_fake[fake_ind], 1, cov=True)
    cov1 = np.sqrt(np.diag(cov1))
    cov2 = np.sqrt(np.diag(cov2))

#    # K-S Statistic
#    Ar = np.array(list(zip(dist_from_cam_real, radius_real)))
#    Af = np.array(list(zip(dist_from_cam_fake, radius_fake)))
#
#    D, p = ks2d.ks2d2s(Ar, Af)
#
#    n = len(Ar)
#    m = len(Af)
#
#    alpha = 2/np.exp(2*D**2/(1/n + 1/m))
#
#    print(D, alpha, n, m)

    x = np.linspace(25, 70, 100)
    y1 = fit1[0]*x + fit1[1]
    y2 = fit2[0]*x + fit2[1]

    ax1.plot(x, y1, c='r', label='$r=(%.3f\pm%.3f)d+(%.1f\pm%.1f)$' % (fit1[0], cov1[0], fit1[1], cov1[1]))
    ax2.plot(x, y2, c='r', label='$r=(%.4f\pm%.4f)d+(%.2f\pm%.2f)$' % (fit2[0], cov2[0], fit2[1], cov2[1]))
    
    ax1.legend()
    ax2.legend()
    
#    s.append(fit1[0])
#    b.append(fit1[1])
    
    if cam == 0:
        ax1.set_title('PICO-60 images')
        ax2.set_title('Artificial images')

fig.text(0.5, 0.07, 'Distance from camera $d$ [cm]', ha='center', size=14)
fig.text(0.045, 0.5, 'Bubble radius $r$ [pixels]', va='center', rotation='vertical', size=14)

fig.text(0.07, 0.8, 'Camera 0', va='center', rotation='vertical', size=12)
fig.text(0.07, 0.6, 'Camera 1', va='center', rotation='vertical', size=12)
fig.text(0.07, 0.4, 'Camera 2', va='center', rotation='vertical', size=12)
fig.text(0.07, 0.2, 'Camera 3', va='center', rotation='vertical', size=12)

fig.text(0.5, 0.92, 'Radius of bubbles in center region of jar', ha='center', size=16)

#plt.show()
#plt.tight_layout()
plt.savefig('radius_at_center.png')

#print(s)
#print(b)

