import numpy as np
import matplotlib.cm as cmx
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib
import matplotlib.colors as colors

#def PolarBarPlot(lv, radD, outfolder, YYYY, DOY, XH, hours, XM, minu, iStep, idxStep):
def PolarBarPlot(lv, solar_altitude, solar_azimuth, fig_title, filename_out, minrad, maxrad, plot_type):
#def PolarBarPlot(lv, filename_start, outfolder, YYYY, DOY, XH, hours, XM, minu, iStep, minrad, maxrad):

    deg2rad = np.pi/180

    fig = figure()
    ax = fig.add_subplot(111, projection = 'polar')

    # Set zero location and clockwise
    ax.set_theta_zero_location("N")
    if plot_type: 
        ax.set_theta_direction('clockwise')
    else:
        ax.set_theta_direction('anticlockwise')
        # ax.set_theta_direction('clockwise')
    skyalt, skyalt_c = np.unique(lv[:, 0], return_counts=True)   # Unique altitudes in lv, i.e. unique altitude for the patches

    # lvSum = np.sum(lv[:,2])

    lv_norm = lv[:, 2]  # Watts per square meter Steradian
    lv_alt = lv[:, 0]
    # lv_azi = lv[:, 1]

    # lvMax = np.around(np.max(lv_norm), decimals=3)

    if plot_type:
        jet = cm = plt.get_cmap('jet')
        cNorm = colors.Normalize(vmin=minrad, vmax=maxrad) # Watts per square meter Steradian
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    else:
        # patch_category = {1.8:'deepskyblue', 2.5:'forestgreen', 3.3:'yellow', 4.5:'peru', 6.0:'grey'}
        patch_category = {1.8:'deepskyblue', 2.5:'forestgreen', 3.3:'yellow', 4.5:'peru'}

    if lv.shape[0] < 160:
        azistart = np.array([0, 4, 2, 5, 8, 0, 10, 0]) # Fredrik/Nils
    else:
        azistart = np.array([0, 0, 4, 4, 2, 2, 5, 5, 8, 8, 0, 0, 10, 10, 0]) # Nils test
    radii = np.ones(np.max(skyalt_c))*1.5
    radii_sub = 0
    for i in range(skyalt_c.__len__()):
        clrs = lv_norm[lv_alt == skyalt[i]]
        # print(clrs)
        if skyalt_c[i] > 1:
            theta_patch = ( np.arange(0, skyalt_c[skyalt == skyalt[i]][0], 1) * (360 / skyalt_c[skyalt == skyalt[i]][0] )) * deg2rad
            width_patch = ( np.ones( skyalt_c[skyalt == skyalt[i]][0] ) * ( 360 * deg2rad ) / skyalt_c[skyalt == skyalt[i]][0] )

            # if plot_type:
            #     patch_order_range = range(skyalt_c[skyalt == skyalt[i]][0])
            # else:
            #     patch_order_range = reversed(range(skyalt_c[skyalt == skyalt[i]][0]))

            patch_order_range = range(skyalt_c[skyalt == skyalt[i]][0])

            for j in patch_order_range:
            # for j in range(skyalt_c[skyalt == skyalt[i]][0]):
                if plot_type:
                    patch_color = scalarMap.to_rgba(clrs[j])
                    # bars = ax.bar(theta_patch[j] + (azistart[i] * deg2rad), radii[j]-radii_sub, width=width_patch[j], bottom=0.0, edgecolor='black', facecolor=patch_color)
                else:
                    patch_color = patch_category[clrs[j]]
                    # bars = ax.bar((np.pi*2 - theta_patch[j]) + (azistart[i] * deg2rad), radii[j]-radii_sub, width=width_patch[j], bottom=0.0, edgecolor='black', facecolor=patch_color)
                # bars = ax.bar(theta_patch[j] + (azistart[i] * deg2rad) - width_patch[j], radii[j]-radii_sub, width=width_patch[j], bottom=0.0, edgecolor='black', facecolor=patch_color)
                bars = ax.bar(theta_patch[j] + (azistart[i] * deg2rad), radii[j]-radii_sub, width=width_patch[j], bottom=0.0, edgecolor='black', facecolor=patch_color)
            if lv.shape[0] < 160:
                radii_sub += 0.2
            else:
                radii_sub += 0.1 # Nils test
        else:
            # Create a circle for the center of the plot
            my_circle=Circle( (0,0), 0.1, edgecolor='black', transform=ax.transData._b)
            if plot_type:
                patch_color = scalarMap.to_rgba(clrs[0])
            else:
                patch_color = patch_category[clrs[0]]
            my_circle.set_facecolor(patch_color)
            ax.add_artist(my_circle)

    # Adding sun's position
    if ((plot_type == 0) & (solar_altitude > 0)):
        sun_position = plt.scatter((solar_azimuth * deg2rad), 1.5 - (1.5 * (solar_altitude/90)), s=(70 * (solar_altitude/90)), color='gold')
        sun_position_star1 = plt.scatter((solar_azimuth * deg2rad), 1.5 - (1.5 * (solar_altitude/90)), s=(400 * (solar_altitude/90)), color='gold', marker='1')
        sun_position_star2 = plt.scatter((solar_azimuth * deg2rad), 1.5 - (1.5 * (solar_altitude/90)), s=(400 * (solar_altitude/90)), color='gold', marker='2')
        # sun_position = plt.scatter(np.pi * 2 - (solar_azimuth * deg2rad), 1.5 - (1.5 * (solar_altitude/90)), s=(70 * (solar_altitude/90)), color='gold')
        # sun_position_star1 = plt.scatter(np.pi * 2 - (solar_azimuth * deg2rad), 1.5 - (1.5 * (solar_altitude/90)), s=(400 * (solar_altitude/90)), color='gold', marker='1')
        # sun_position_star2 = plt.scatter(np.pi * 2 - (solar_azimuth * deg2rad), 1.5 - (1.5 * (solar_altitude/90)), s=(400 * (solar_altitude/90)), color='gold', marker='2')
        ax.add_artist(sun_position_star1)
        ax.add_artist(sun_position_star2)
        ax.add_artist(sun_position)

    # Remove grids and tick labels
    ax.set_yticklabels([])
    ax.grid(False)
    ax.spines['polar'].set_visible(False)

    if plot_type:
        # Adding colorbar
        sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=minrad, vmax=maxrad)) # Watts per square meter Steradian
        sm._A = []
        cbaxes = fig.add_axes([0.87, 0.1, 0.03, 0.8])
        cb = plt.colorbar(sm, cax = cbaxes)
        cb.ax.set_title(r'$W/m^2$ $sr^{-1}$', fontsize=12, fontweight='bold') # Watts per square meter Steradian
    else:
        # legend_colors = {'Shaded building wall':'peru', 'Sunlit building wall':'yellow', 'Sky':'deepskyblue', 'Trees':'forestgreen', 'Roof':'black'}
        # legend_colors = {'Building wall':'peru', 'Building roof':'grey', 'Trees':'forestgreen', 'Sky':'deepskyblue'}
        legend_colors = {'Building wall':'peru', 'Trees':'forestgreen', 'Sky':'deepskyblue'}
        legend_labels = list(legend_colors.keys())
        legend_handles = [plt.Rectangle((0,0),1,1, color=legend_colors[legend_label]) for legend_label in legend_labels]
        if solar_altitude > 0:
            legend_handles.extend([(sun_position, sun_position_star1, sun_position_star2)])
            legend_labels.extend(['Sun'])
        plt.legend(legend_handles, legend_labels,
                    title='Categories',
                    loc='lower left',
                    bbox_to_anchor=(0.95, 0, 0.5, 1),
                    fontsize=7,
                    frameon=False,
                    title_fontsize=7,
                    markerscale=0.7)

    ax.set_title(fig_title, weight='bold', fontsize=10)

    fig.savefig(filename_out, dpi=300)
    plt.close('all')