import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm
import csv

from gmt import write_gmt



def create_legend(color_dict, cytoscape_file_tag):
    # make a legend
    fig, ax = plt.subplots()

    #create a grid for positioning the rectangles
    grid = np.array([np.repeat(0.3, 2 * len(color_dict.keys())), np.flip(np.linspace(0.1,0.9,2*len(color_dict.keys())))]).T

    
    # rectangle dimension
    rec_dim = 0.1 / (len(color_dict.keys()) / 3)
    label_font_size = 12 / (len(color_dict.keys()) / 3)
        
    ax.annotate('Contribution', np.array(grid[0]) + np.array([0.15,0]), color='black', weight='bold', 
                fontsize=label_font_size, ha='center', va='center')

    # create rectangles
    for i, label in enumerate(list(color_dict.keys())):
        # define rectangle sizes
        rect = mpatches.Rectangle(grid[i + 1], rec_dim, rec_dim, color = color_dict[label])
        ax.add_patch(rect)
        
        # annotate
        ax.add_artist(rect)
        rx, ry = rect.get_xy()
        cx = rx + rect.get_width()/2
        cy = ry + rect.get_height()/2
        ax.annotate(label, (cx + 0.15, cy), color='black', weight='bold', 
                    fontsize=label_font_size, ha='center', va='center')



    plt.axis('off')
    outfile = cytoscape_file_tag + 'legend.pdf'
    plt.savefig(outfile)

def prepareCytoscape(terms, gmt, cytoscape_file_tag, col_significance = None, color_palette = None, custom_colors = None, color_integrated_only="#FFFFF0"):

    if type(col_significance) != None:
        # Obtain the name of each omics dataset and incorporate a 'combined' contribution
        tests = col_significance.columns.to_numpy(dtype='<U64')[2:]
        tests = list(map(lambda x: x[6:], tests))
        tests.append("combined")

        # Create a matrix of ones and zeros, where columns are omics datasets + 'combined' and rows are enriched pathways
        evidence_cols = pd.DataFrame()
        for i in range(len(tests)):
            evidence_cols.loc[:,tests[i]] = list(map(lambda x: int(tests[i] in x), col_significance['evidence']))
        #evidence_cols.loc[:,'term_id'] = col_significance['term_id'].to_numpy()
        evidence_cols.insert(loc=0, column='term_id', value=col_significance['term_id'].to_numpy())
        
        mapped_colors = []
        if color_palette == None and custom_colors == None:
            # Create default color palette
            cmap = mpl.cm.get_cmap('Set1', len(tests)-1)
            mapped_colors = list(map(lambda x: mpl.colors.rgb2hex(cmap(x)),range(len(tests)-1)))
            
        elif custom_colors != None:
            mapped_colors = custom_colors   
            
        # Create user-defined color pallette
        else:
            cmap = mpl.cm.get_cmap(color_palette, len(tests)-1)
            mapped_colors = list(map(lambda x: mpl.colors.rgb2hex(cmap(x)),range(len(tests)-1)))
            
        mapped_colors.append(color_integrated_only)

        # create a dictionary for colors
        color_dict = {}
        for i, color in enumerate(mapped_colors):
            color_dict[tests[i]] = color

        # cytoscape instructions to create pie charts
        instruct_str = f'''piechart: attributelist="{",".join(tests)}" colorlist="{",".join(mapped_colors)}" showlabels=FALSE'''
        evidence_cols.loc[:,'instruct'] = instruct_str
        
        # Writing the Files
        outfile = cytoscape_file_tag + "subgroups.txt"
        evidence_cols.to_csv(outfile, sep='\t', index=False, quoting=csv.QUOTE_NONE)

        create_legend(color_dict, cytoscape_file_tag)

    # if there is only one dataset, write results to files
    outfile = cytoscape_file_tag + "pathways.txt"
    terms.to_csv(outfile, sep='\t', index=False)

    write_gmt(gmt, (cytoscape_file_tag + "pathways.gmt"))
