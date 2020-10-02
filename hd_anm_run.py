import argparse
import numpy
import urllib.request

from packman import molecule
from packman.anm import hdANM
from packman.constants import amino_acid_molecular_weight
'''
##################################################################################################
#                                          Interface                                             #
##################################################################################################
'''

def IO():
    """User interface for the user to provide the parameters
    Todo:
        * None
    
    Returns:
        Namespace: Various arguments in various formats
    """
    parser=argparse.ArgumentParser(description='Rigid Domain ANM (RD-ANM). (https://github.com/Pranavkhade/PACKMAN)')
    parser.add_argument('-pdbid','--pdbid', metavar='PDB_ID', type=str, help='If provided, the PBD with this ID will be downloaded and saved to FILENAME.')
    parser.add_argument('filename', metavar='FILENAME', help='Path and filename of the PDB file.')
    parser.add_argument('hngfile', metavar='HNG', help='Path and filename of the corresponding HNG file.')

    parser.add_argument("--chain", help='Enter The Chain ID')
    parser.add_argument("--dr", type=float, default=15, help='Distance cutoff for the ANM.')
    parser.add_argument("--power", type=float, default=0, help='Power of the distance in non-parametric ANM.')
    parser.add_argument("--mass", default='residue', help='Mass of the residue; unit or molecular weight')


    parser.add_argument("--scale", type=int, default=2, help='movie scale')
    parser.add_argument("--frames", type=int, default=10, help='number of frames')
    parser.add_argument("--modes", type=int, default=10, help='how many modes')

    # web server parameters
    web_server_group = parser.add_argument_group('Web server parameters', 'Used by the web form')
    web_server_group.add_argument('--callbackurl', type=str, help='Optional callback url if this script was called from Drupal.')
    web_server_group.add_argument('--nodeid', type=int, help='Optional node id if this script was called from Drupal.')

    args=parser.parse_args()
    return args


'''
##################################################################################################
#                                              Main                                              #
##################################################################################################
'''


def main():
    """
    """
    ARGS = IO()


    try:
        #Load the Hinge Information
        filename = ARGS.filename
        chain    = ARGS.chain
        dr       = float(ARGS.dr)
        power    = float(ARGS.power)

        ##@ResearchIT PAGE/ SECTION 1 (Title: Input File)

        #@ResearchIT, this is the part where we give user an option to either download or upload the file so next line will be optional if user selects to upload the file.
        #mol=molecule.download_structure(ARGS.pdbid,filename)

        mol=molecule.load_structure(filename)

        #@ResearchIT Chain, if specified, if not, All
        #Message on the website: This server allows C-alpha atoms only because of the computational constraints; however, the all-atom model can run using PACKMAN - hdANM API on the user's computer. Please refer to the following page.
        #https://py-packman.readthedocs.io/en/latest/tutorials/hdANM.html#tutorials-hdanm
        if(chain is not None):
            calpha=[i for i in mol[0][chain].get_calpha() if i is not None]
        else:
            calpha=[i for i in mol[0].get_calpha()  if i is not None]


        ##@ResearchIT PAGE/ SECTION 2
        #Arguments
        #dr= cutoff value (float)
        #power = Power of Distance (float)
        #Mass = Unit or Molecular Weight of the Residue (only two options so dropdown miight work)

        #Message of the website: hdANM also requires the user to provide the .hng file containing the information about which parts of the protein are hinges and which parts are rigid domains. This information can be obtained from the PACKMAN’s hinge prediction algorithm. Please read the following for more details: https://packman.bb.iastate.edu/
        #Hinge Information file format (.hng) Example:
        #Code block:
        #1exr.pdb_A  D1  1:70
        #1exr.pdb_A  H1  71:90
        #1exr.pdb_A  D2  91:148
        #The first column in the .hng file is Filename_ChainID, the second column is Domain/Hinge ID, and the third column is the residues in the particular domain/hinge. The .hng file is tab-separated. (Tabs separate columns)


        #@ResearchIT, is it possible to give an optino to add hinge and domains from the interface line by line (optional to the file upload)
        #eg... user given dropdown box saying 'add new hinge'/ 'add new domain'. If user selects 'add new hinge', they will get two text boxes asking for filename_chainID (first column) and residue range (last column) whereas hinge id will be calculated automatically H1,H2.... so on. and same for the domain.
        

        Model=hdANM(calpha,dr=dr,power=power,hng_file=ARGS.hngfile)
        Model.calculate_hessian(mass_type=ARGS.mass)
        Model.calculate_decomposition()

        
        numpy.savetxt("eigenvalues.csv", Model.get_eigenvalues(), delimiter=",")
        numpy.savetxt("eigenvectors.csv", Model.get_eigenvectors(), delimiter=",")

        

        ##@ResearchIT PAGE/ SECTION 3
        #For more details on the curvilinear projection (movie) please read: 
        #parameters here
        #ARGS.modes: How many modes do you want to get a curvilinear projection (movie)? (Note: First non-rigid mode is 6.pdb; Second non-rigid mode is 67.pdb and so on)
        #Scale: Scale for the movie
        #Frames: Number of frames you want in the movies.
        #Message on Website: Eigenvalues and Eigenvector indecies start from 0 (includes rigid body motion)


        for i in range(6,6+ARGS.modes,1):
            Model.calculate_movie(i,scale=ARGS.scale,n=ARGS.frames)
        
        if ARGS.nodeid is not None:
            urllib.request.urlopen(ARGS.callbackurl + '/' + str(ARGS.nodeid) + "/0")

    except Exception:
        if ARGS.nodeid is not None:
            urllib.request.urlopen(ARGS.callbackurl + '/' + str(ARGS.nodeid) + "/1")


    return True


'''
##################################################################################################
#                                               IO                                               #
##################################################################################################
'''

if(__name__=='__main__'):
    #@ResearchIt: how to run : python.exe .\hd_anm_run.py .\1prw.pdb .\output.hng.txt
    main()
