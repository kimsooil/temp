
from sklearn.externals import *

#import pickle
import joblib

from rdkit import Chem # chemoinformatics package
from rdkit.Chem import Draw

import deepchem as dc # ML package 
from deepchem.feat import CircularFingerprint # featurizer

import os, time, datetime
import numpy
from flask import Flask, render_template, request, jsonify, make_response

from pymatgen import Composition

app = Flask(__name__)

Hartree=27.2114 # 1 Ha = 27.2114 eV
#path = "/home/kimsooil/matrics/" # to be changed!!
path = "./" # to be changed!!

def getInchiFromSmiles_rdkit(smiles='C1(/C=C)=C\C=C/C=C1'):
    return Chem.MolToInchi(Chem.MolFromSmiles(smiles))

@app.route("/convert_smiles_to_inchi", methods=['GET'])
def getInchi():
    smiles = request.args.get('smiles')
    inchi=""
    inchi=getInchiFromSmiles_rdkit(smiles)
    return jsonify({"smiles":smiles, "inchi":inchi})
    
#def predict_gap(smiles='C1(/C=C)=C\C=C/C=C1'):
def predict_gap(smiles='C1(/C=C)=C\C=C/C=C1'):
    if not os.path.exists('/tmp/tmpujbox9yr'):
        os.makedirs('/tmp/tmpujbox9yr')

    #filename = path + "qm9db_gap-predictor_RF-100estimators-0.5maxfeatures.mod.joblib" # Hartree?
    filename = path + "qm9+cep_gap-predictor_RF-50estimators-0.5maxfeatures.mod0.joblib" # ev? send as _Ha so it converted to eV later

    model = joblib.load(filename)
    # define featurizer (Extended Connectivity Finger-Print)
    ecfp = dc.feat.CircularFingerprint(size=1024, radius=2)
    smiles=smiles.upper()
    test_smiles =[smiles]

    # convert SMILE to 2D molecule
    test_mol = [Chem.MolFromSmiles(s) for s in test_smiles]
    print(test_mol)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S.%f")
    svg=Chem.Draw.MolToFile(test_mol[0], 'static/organic-ml-predict'+timestamp+'.png')
    natoms=0
    inchi=""
    if test_mol[0]:
        #natoms=Chem.AddHs(test_mol[0]).GetNumAtoms()
        inchi=Chem.MolToInchi(test_mol[0])
        formula=inchi.split('/')
        formula=formula[1]
        natoms=Composition(formula).num_atoms
        print(inchi, natoms)
    # featurize molecule
    if test_mol[0] is not None:
        test_ecfp = ecfp.featurize(test_mol)
        # predict gap (in Hartree)
        predicted_gap = model.predict_on_batch(test_ecfp)
        return natoms, predicted_gap[0], inchi, timestamp
    else:
        return -1

#@app.route("/homo", methods=['GET'])
def predict_homo(smiles):

    smiles = request.args.get('smiles')

    if not os.path.exists('/tmp/tmpujbox9yr'):
        os.makedirs('/tmp/tmpujbox9yr')

    # load model from disk
    
    filename = path + "qm9_homo-predictor_RF-100estimators-0.5maxfeatures.mod.joblib" # ev? send as _Ha so it converted to eV later

    model = joblib.load(filename)
    # define featurizer (Extended Connectivity Finger-Print)
    ecfp = dc.feat.CircularFingerprint(size=1024, radius=2)
    smiles=smiles.upper()
    test_smiles =[smiles]

    # convert SMILE to 2D molecule
    test_mol = [Chem.MolFromSmiles(s) for s in test_smiles]

    # featurize molecule
    predicted=-1000
    if test_mol[0] is not None:
        test_ecfp = ecfp.featurize(test_mol)
        # predict gap (in Hartree)
        predicted = model.predict_on_batch(test_ecfp)
        predicted=predicted[0]

    #return jsonify({"homo":predicted, "smiles":smiles})
    return predicted

#@app.route("/lumo", methods=['GET'])
def predict_lumo(smiles):

    smiles = request.args.get('smiles')

    if not os.path.exists('/tmp/tmpujbox9yr'):
        os.makedirs('/tmp/tmpujbox9yr')

    # load model from disk
    
    filename = path + "qm9_lumo-predictor_RF-100estimators-0.5maxfeatures.mod.joblib" # ev? send as _Ha so it converted to eV later

    model = joblib.load(filename)
    # define featurizer (Extended Connectivity Finger-Print)
    ecfp = dc.feat.CircularFingerprint(size=1024, radius=2)
    smiles=smiles.upper()
    test_smiles =[smiles]

    # convert SMILE to 2D molecule
    test_mol = [Chem.MolFromSmiles(s) for s in test_smiles]

    # featurize molecule
    predicted=-1000
    if test_mol[0] is not None:
        test_ecfp = ecfp.featurize(test_mol)
        # predict gap (in Hartree)
        predicted = model.predict_on_batch(test_ecfp)
        predicted=predicted[0]

    #return jsonify({"lumo":predicted, "smiles":smiles})
    return predicted

#@app.route("/zpve", methods=['GET'])
def predict_zpve(smiles):

    smiles = request.args.get('smiles')

    if not os.path.exists('/tmp/tmpujbox9yr'):
        os.makedirs('/tmp/tmpujbox9yr')

    # load model from disk
    
    filename = path + "qm9_zpve-predictor_RF-100estimators-0.5maxfeatures.mod.joblib" # ev? send as _Ha so it converted to eV later

    model = joblib.load(filename)
    # define featurizer (Extended Connectivity Finger-Print)
    ecfp = dc.feat.CircularFingerprint(size=1024, radius=2)
    smiles=smiles.upper()
    test_smiles =[smiles]

    # convert SMILE to 2D molecule
    test_mol = [Chem.MolFromSmiles(s) for s in test_smiles]

    # featurize molecule
    predicted=-1000
    if test_mol[0] is not None:
        test_ecfp = ecfp.featurize(test_mol)
        # predict gap (in Hartree)
        predicted = model.predict_on_batch(test_ecfp)
        predicted=predicted[0]

    #return jsonify({"zpve":predicted, "smiles":smiles})
    return predicted

#@app.route("/gap", methods=['GET'])
def predict_gap_only(smiles):

    smiles = request.args.get('smiles')

    if not os.path.exists('/tmp/tmpujbox9yr'):
        os.makedirs('/tmp/tmpujbox9yr')

    # load model from disk
    
    filename = path + "qm9+cep_gap-predictor_RF-50estimators-0.5maxfeatures.mod0.joblib" # ev? send as _Ha so it converted to eV later

    model = joblib.load(filename)
    # define featurizer (Extended Connectivity Finger-Print)
    ecfp = dc.feat.CircularFingerprint(size=1024, radius=2)
    smiles=smiles.upper()
    test_smiles =[smiles]

    # convert SMILE to 2D molecule
    test_mol = [Chem.MolFromSmiles(s) for s in test_smiles]

    # featurize molecule
    predicted=0
    if test_mol[0] is not None:
        test_ecfp = ecfp.featurize(test_mol)
        # predict gap (in Hartree)
        predicted = model.predict_on_batch(test_ecfp)
        predicted=predicted[0]

    #return jsonify({"gap":predicted, "smiles":smiles})
    return predicted

@app.route("/", methods=['GET'])
def hello_world():
    start=time.time()
    smiles = request.args.get('smiles')
    natoms, predicted_gap, inchi, timestamp=predict_gap(smiles)
    #predicted_gap_Ha = predicted_gap / Hartree # comment out when unit is eV
    p_homo=predict_homo(smiles)
    #p_gap=predict_gap_only(smiles)
    #p_lumo=predict_lumo(smiles)
    p_zpve=predict_zpve(smiles)
    return jsonify({"msg":"Hello World", "smiles":smiles.upper(), "natoms":natoms, "inchi":inchi, "timestamp":timestamp,
        "unit":"Hartree", "predicted_gap":predicted_gap,"homo":p_homo,"lumo":p_homo+predicted_gap, "gap":predicted_gap, "zpve":p_zpve,
        "time_consumed":time.time()-start })

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=7700, debug=True, threaded=True)
