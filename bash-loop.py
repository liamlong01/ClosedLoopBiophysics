
import sys
import shutil
import os
import run
import pickle

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.collections import LineCollection

dllfile = 'nrnmech.dll'

home = "D:/OneDrive - University of Toronto/cns/sims"
results = "E:/results"
start_T = 0
end_T = 3000
dt = 0.025

def plotstuff(cell, electrode):
    '''plotting'''
    fig = plt.figure(dpi=160)

    ax1 = fig.add_axes([0.05, 0.1, 0.55, 0.9], frameon=False)
    cax = fig.add_axes([0.05, 0.115, 0.55, 0.015])

    ax1.plot(electrode.y, electrode.z, '.', marker='o', markersize=1, color='k',
             zorder=0)

    # normalize to min peak
    LFPmin = electrode.LFP.min(axis=1)
    LFPnorm = -(electrode.LFP.T / LFPmin).T

    i = 0
    zips = []
    for x in LFPnorm:
        zips.append(list(zip(cell.tvec /1000 + electrode.x[i] + 2,
                             x * 12 + electrode.z[i])))
        i += 1

    line_segments = LineCollection(zips,
                                   linewidths=(1),
                                   linestyles='solid',
                                   cmap='nipy_spectral',
                                   zorder=1,
                                   rasterized=False)
    line_segments.set_array(np.log10(-LFPmin))
    ax1.add_collection(line_segments)

    axcb = fig.colorbar(line_segments, cax=cax, orientation='horizontal')
    axcb.outline.set_visible(False)
    xticklabels = np.array([-0.1, -0.05, -0.02, -0.01, -0.005, -0.002])
    xticks = np.log10(-xticklabels)
    axcb.set_ticks(xticks)
    axcb.set_ticklabels(np.round(-10 ** xticks, decimals=3))
    axcb.set_label('spike amplitude (mV)', va='center')

    ax1.plot([22, 38], [100, 100], color='k', lw=1)
    ax1.text(22, 102, '10 ms')

    ax1.plot([60, 80], [100, 100], color='k', lw=1)
    ax1.text(60, 102, '20 $\mu$m')

    ax1.set_xticks([])
    ax1.set_yticks([])

    axis = ax1.axis(ax1.axis('equal'))
    ax1.set_xlim(axis[0] * 1.02, axis[1] * 1.02)

    # plot morphology
    zips = []
    for y, z in cell.get_pt3d_polygons(projection=('y','z')):
        zips.append(list(zip(y, z)))
    from matplotlib.collections import PolyCollection
    polycol = PolyCollection(zips, edgecolors='none',
                             facecolors='gray', zorder=-1, rasterized=False)
    ax1.add_collection(polycol)

    ax1.text(-0.05, 0.95, 'a',
             horizontalalignment='center',
             verticalalignment='center',
             fontsize=16, fontweight='demibold',
             transform=ax1.transAxes)

    # plot extracellular spike in detail
    ind = np.where(electrode.LFP == electrode.LFP.min())[0][0]
    timeind = (cell.tvec >= 0) & (cell.tvec <= 10)
    xticks = np.arange(10)
    xticklabels = xticks
    LFPtrace = electrode.LFP[ind,]
    vline0 = cell.tvec[cell.somav == cell.somav.max()]
    vline1 = cell.tvec[LFPtrace == LFPtrace.min()]
    vline2 = cell.tvec[LFPtrace == LFPtrace.max()]

    # plot asterix to link trace in (a) and (c)
    ax1.plot(electrode.x[ind], electrode.z[ind], '*', markersize=5,
             markeredgecolor='none', markerfacecolor='k')

    ax2 = fig.add_axes([0.75, 0.6, 0.2, 0.35], frameon=True)
    ax2.plot(cell.tvec[timeind], cell.somav[timeind], lw=1, color='k', clip_on=False)

    ax2.vlines(vline0, cell.somav.min(), cell.somav.max(), 'k', 'dashed', lw=0.25)
    ax2.vlines(vline1, cell.somav.min(), cell.somav.max(), 'k', 'dashdot', lw=0.25)
    ax2.vlines(vline2, cell.somav.min(), cell.somav.max(), 'k', 'dotted', lw=0.25)

    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticks)
    ax2.axis(ax2.axis('tight'))
    ax2.set_ylabel(r'$V_\mathrm{soma}(t)$ (mV)')

    for loc, spine in ax2.spines.items():
        if loc in ['right', 'top']:
            spine.set_color('none')
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')

    ax2.set_title('somatic potential', va='center')

    ax2.text(-0.3, 1.0, 'b',
             horizontalalignment='center',
             verticalalignment='center',
             fontsize=16, fontweight='demibold',
             transform=ax2.transAxes)

    ax3 = fig.add_axes([0.75, 0.1, 0.2, 0.35], frameon=True)
    ax3.plot(cell.tvec[timeind], LFPtrace[timeind], lw=1, color='k', clip_on=False)
    ax3.plot(0.5, 0, '*', markersize=5, markeredgecolor='none', markerfacecolor='k')

    ax3.vlines(vline0, LFPtrace.min(), LFPtrace.max(), 'k', 'dashed', lw=0.25)
    ax3.vlines(vline1, LFPtrace.min(), LFPtrace.max(), 'k', 'dashdot', lw=0.25)
    ax3.vlines(vline2, LFPtrace.min(), LFPtrace.max(), 'k', 'dotted', lw=0.25)

    ax3.set_xticks(xticks)
    ax3.set_xticklabels(xticks)
    ax3.axis(ax3.axis('tight'))

    for loc, spine in ax3.spines.items():
        if loc in ['right', 'top']:
            spine.set_color('none')
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')

    ax3.set_xlabel(r'$t$ (ms)', va='center')
    ax3.set_ylabel(r'$\Phi(\mathbf{r},t)$ (mV)')

    ax3.set_title('extracellular spike', va='center')

    ax3.text(-0.3, 1.0, 'c',
             horizontalalignment='center',
             verticalalignment='center',
             fontsize=16, fontweight='demibold',
             transform=ax3.transAxes)

    return fig

def create_stimuli(sec, step_number):
    """Create the stimuli"""

    print('Attaching stimulus electrodes')

    stimuli = []
    step_amp = [0] * 3

    with open('current_amps.dat', 'r') as current_amps_file:
        first_line = current_amps_file.read().split('\n')[0].strip()
        hyp_amp, step_amp[0], step_amp[1], step_amp[2] = first_line.split(' ')

    iclamp = neuron.h.IClamp(0.5, sec=sec)
    iclamp.delay = 700
    iclamp.dur = 2000
    iclamp.amp = float(step_amp[step_number - 1])
    print('Setting up step current clamp: '
          'amp=%f nA, delay=%f ms, duration=%f ms' %
          (iclamp.amp, iclamp.delay, iclamp.dur))

    stimuli.append(iclamp)

    hyp_iclamp = neuron.h.IClamp(0.5, sec=sec)
    hyp_iclamp.delay = 0
    hyp_iclamp.dur = 3000
    hyp_iclamp.amp = float(hyp_amp)
    print('Setting up hypamp current clamp: '
          'amp=%f nA, delay=%f ms, duration=%f ms' %
          (hyp_iclamp.amp, hyp_iclamp.delay, hyp_iclamp.dur))

    stimuli.append(hyp_iclamp)

    return stimuli


def getTemplate(directory = ''):
    f = open(directory + 'template.hoc')
    for line in f:
        if 'begintemplate' in line:
            templatename = line.split(' ')[1][:-1]
            #print(templatename)
            return templatename


if __name__ == '__main__':

    neuronfolder = 'neurons/' + sys.argv[1] +'/'

# for retry in range(100):
#     try:
#         os.remove('nrnmech.dll')
#         break
#     except PermissionError as e:
#         print(e)
# else:
#     print("could not remove nrn mech file (e.g.) permissions")
#     exit(1)

# if 'nrnmech.dll' not in os.listdir():
#     if 'nrnmech.dll' in os.listdir(neuronfolder + 'mechanisms'):
#         shutil.copy(neuronfolder + 'mechanisms/nrnmech.dll', 'nrnmech.dll')
#     else: #TODO compile neuron mechanisms if file is not found
#         print("Compiled mechanisms not find. TODO try to compile if this happens")
#         exit(1)

    # we need to delay these imports until we get the correct mechanism file
	#TODO move these imports? doing mechanism loading in bash now

    import LFPy
    import neuron

    #cwd = os.getcwd()
    #os.chdir(neuronfolder)



    # delete old sections from NEURON namespace
    LFPy.cell.neuron.h("forall delete_section()")



    neuron.h.load_file("stdrun.hoc")
    neuron.h.load_file("import3d.hoc")

    print('Loading constants')
    neuron.h.load_file('constants.hoc')

    neuron.h.load_file(1, "morphology.hoc")
    neuron.h.load_file(1, "biophysics.hoc")

    print('Loading constants')
    neuron.h.load_file(1,'template.hoc')
    neuron.h.load_file('synapses/synapses.hoc')

    print("running.....")

    morphologyfile = 'morphology/' + os.listdir('morphology')[0]  # glob('morphology\\*')[0]
    
    templatename = getTemplate(directory = '')
    


    cell = LFPy.TemplateCell(morphology= morphologyfile,
                             templatefile='template.hoc',
                             templatename=templatename,
                             templateargs=0,
                             tstop=end_T,
                             tstart=start_T,
                             dt=dt,
                             v_init= -70,
                             pt3d=True,
                             delete_sections=True,
                             verbose=True)

    for sec in cell.somalist:
        somasec = sec

    stimuli = create_stimuli(somasec, 1)



    res = 10



    xmin = min(cell.xend)


    ymin = int(min(cell.yend)/res) - 5
    ymax = int(max(cell.yend)/res) + 5

    zmin = int(min(cell.zend)/res) - 5
    zmax = int(max(cell.zend)/res) + 5
    """
    # Generate the grid in xz-plane over which we calculate local field potentials
    X, Y, Z = np.mgrid[1:2, ymin:ymax:1,   zmin:zmax:1] * res
    X = X - xmin
    # define parameters for extracellular recording electrode, using optional method
    electrodeParameters = {
        'sigma': 0.3,  # extracellular conductivity
        'x': X.flatten(),  # x,y,z-coordinates of contacts
        'y': Y.flatten(),
        'z': Z.flatten(),
        'method': 'soma_as_point',  # sphere source soma segment
        'N': np.array([[1, 0, 0]] * X.size),  # surface normals
        'r': 6,  # contact site radius
        'n': 20,  # datapoints for averaging
    }

    print("creating electrode...")

    print("z-direction from : %d %d" % (zmin,zmax))
    print("y-direction from : %d %d" % (ymin,ymax))

    print("total electrodes: %d" % len(X.flatten()))
    # create extracellular electrode object for LFPs on grid
    electrode = LFPy.RecExtElectrode(**electrodeParameters)

    """
	
    print( "simulating....")

    cell.simulate(rec_imem=True)

    print("done simulating...")

    os.chdir(home  )
    if not os.path.exists('results/' ):
        os.makedirs('results/')
    os.chdir(results)


    print("pickling cell..")
    neurontype = sys.argv[1].replace('/','')
    try:
        cell.cellpickler('cellvoltages_'+ neurontype)
    except Exception as e:
        print("Error pickling cell:")
        print(e)
		
    """
    print("pickling electrode..")
    try:
        f = open('plenty_electrodelfps_' + neurontype, 'wb')
        pickle.dump(electrode, f)

    except Exception as e:
        print("Error pickling electrode:")
        print(e)

    f.close()
    """
	
	
    os.chdir(home)
    exit(0)
  #  os.chdir(home)

   # plt.plot(cell.tvec, cell.somav)

    #fig = plt.figure()

    #plt.plot(cell.tvec, electrode.LFP[100])



    #plt.show()












