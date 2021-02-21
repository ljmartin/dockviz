import streamlit as st

import pandas as pd
import numpy as np
import os

from rdkit import Chem, RDConfig
from rdkit.Chem import Draw, ChemicalFeatures, Descriptors, Crippen

import py3Dmol


st.set_page_config(
    layout='wide',
    )

##Setup:

def isRingAromatic(mol, atomRing):
        for id in atomRing:
            if not mol.GetAtomWithIdx(id).GetIsAromatic():
                return False
        return True

def calcAromCOM(pdb, ring):
    pos = pdb.GetConformer(0).GetPositions()
    return pos[np.array(ring)].mean(0)

col = {'Donor':'green', 'Acceptor':'purple'}

def main():

    #print out some explanation stuff in the sidebar:
    #st.sidebar.title("Binding site viz with py3dmol")
    #st.sidebar.write("instructions")


    #and some intro text in the main frame:
    st.title('Docking viz for the D4 receptor using py3dmol')
    st.write("The co-crystallized ligand is present in blue. Input a docking rank (0-4) to see the best-scoring docking conformation, in orange.")
    st.write("Shift-click to zoom, and rotate with primary click")
    labels = st.radio("Annotate with pharmacophore labels", ('Yes', 'No'))
    
    selected_rank = st.text_input("Input Rank (i.e. `4`, ranks are zero-indexed)", "")
    if len(selected_rank)>=1:
            pdb = Chem.MolFromPDBFile('./files/'+selected_rank+'.pdb')
            molwt = Descriptors.MolWt(pdb)
            clogp = Crippen.MolLogP(pdb)
            formalcharge = sum([i.GetFormalCharge() for i in pdb.GetAtoms()])
            st.write(f"**Rank**: {selected_rank}, **MW**: {molwt}, **cLogP**: {clogp}, **Formal Charge**: {formalcharge}")
            Draw.MolToFile(Chem.MolFromSmiles(Chem.MolToSmiles(pdb)), 'mol.png')
            st.image('mol.png')

        
    ###now the app:
    #camera:
    viewer = py3Dmol.view(width=1200, height=800)
    viewer.setBackgroundColor(0x000000)
    viewer.setCameraParameters({ 'fov': '50', 'z':150 });
    
    
    #molecules:
    viewer.addModel(open('./files/5WIU-CAC.pdb', 'r').read(), 'pdb')
    viewer.zoomTo()
    viewer.addModel(open('./files/5WIU-receptor.pdb', 'r').read(), 'pdb')
    if len(selected_rank) >= 1:
        viewer.addModel(open('./files/'+selected_rank+'.pdb', 'r').read(), 'pdb')
        
        
        ri = pdb.GetRingInfo()
        for ring in ri.AtomRings():
            if isRingAromatic(pdb, ring):
                com = calcAromCOM(pdb, ring)
                viewer.addSphere({'radius': 1.25, 
                          'center': {'x':com[0], 'y':com[1], 'z':com[2]}, 
                          'wireframe':'True', 
                            'color':'brown'})
                if labels=='Yes':
                        viewer.addLabel('Aromatic', {'position': {'x':com[0], 'y':com[1], 'z':com[2]}, 'backgroundColor': "brown", 'backgroundOpacity': 0.8})

        #features:
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
        feats = factory.GetFeaturesForMol(pdb)
        for f in [f for f in feats if (f.GetFamily()=='Donor' or f.GetFamily()=='Acceptor')]:
                p = np.array(f.GetPos())
                viewer.addSphere({'radius': 0.7, 
                          'center': {'x':p[0], 'y':p[1], 'z':p[2]}, 
                          'wireframe':'False', 
                          'color':'black',
                                  'opacity':0.8})
                offset = p + {'Donor':0.1, 'Acceptor':-0.1}[f.GetFamily()]
                if labels=='Yes':
                        viewer.addLabel(f.GetFamily(), {'position': {'x':offset[0], 'y':offset[1], 'z':offset[2]}, 'backgroundColor': col[f.GetFamily()], 'backgroundOpacity': 0.8});
        positions = pdb.GetConformer(0).GetPositions()
        for p in [positions[count] for count, i in enumerate(pdb.GetAtoms()) if i.GetFormalCharge()!=0]:
                viewer.addSphere({'radius': 0.8, 
                          'center': {'x':p[0], 'y':p[1], 'z':p[2]}, 
                          #'wireframe':'True', 
                          'color':'#00FF00','opacity':0.8})
                if labels=='Yes':
                        viewer.addLabel("Charge", {'position': {'x':p[0], 'y':p[1], 'z':p[2]}, 'backgroundColor': "pink", 'backgroundOpacity': 0.8})
        
        
    ##selections and stylings:
    #binding residues:
    bresis=[87, 88, 90, 91, 94, 101, 111, 112,115, 116, 117,
            119, 120, 123, 185, 186, 187, 192, 193, 194, 196,
            197, 200, 407, 410, 411, 414, 434, 438,]
    #CCR5
    prot = {'resn': ["AQD", "UNL"], 'invert': 1}
    viewer.setStyle(prot, {'cartoon':{'colorscheme':'cyanCarbon'}}) 
    viewer.addSurface(py3Dmol.SES,{'opacity':0.9, 'color': 'white'}, prot)
    
    viewer.setStyle({'resi':bresis}, {'stick': {'colorscheme': 'blueCarbon'}})
    
    #maraviroc:
    viewer.setStyle({'resn':['AQD']}, {'stick': {'colorscheme': 'blueCarbon'}})
    
    #docked ligands:
    viewer.setStyle({'resn':['UNL']}, {'stick': {'colorscheme': 'orangeCarbon'}})
    #viewer.addSurface(py3Dmol.SES, {'opacity': 0.7}, {'resn':['UNL']})


    #view it;
    viewer.render()
    

    t =viewer.js()
    f = open('viz.html', 'w')
    f.write(t.startjs)
    f.write(t.endjs)
    f.close()

    st.components.v1.html(open('viz.html', 'r').read(), width=1200, height=800)


    st.write("""## Credits:""")
    st.write("This uses `py3dmol` for 3d visualization: https://3dmol.csb.pitt.edu/ , Rego & Koes: https://doi.org/10.1093/bioinformatics/btu829")
    st.write("Molecular properties are calculated with the `RDKit`: https://www.rdkit.org/")
    st.write("Docking used `smina`: https://sourceforge.net/projects/smina/")
    st.write("Docking pipeline was described by `Esben Jannik Bjerrum` at https://www.cheminformania.com/ligand-docking-with-smina/")
    st.write("Ligands are from a publicly available dataset from `Lyu et al`: https://doi.org/10.1038/s41586-019-0917-9")
    st.write("and these pieces were put together by lewis martin: https://ljmartin.github.io/")
    
if __name__=="__main__":
    main()
