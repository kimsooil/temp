import os, time, datetime, pickle
import numpy as np
from pymatgen import Composition, Element, MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, PourbaixEntry
import pymongo
import itertools

def printtimestamp(label=""):
    ts = time.time()
    print(label+':'+datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

def getAllCombinations(inlist):
    outlist = []
    for i in range(1, len(inlist)+1):
        outlist.extend(itertools.combinations(inlist, i))
    #outlist = [x[0] if len(x) == 1 else list(x) for x in outlist]
    outlist = [list(x) for x in outlist]
    return outlist

def getPourbaixEntryfromMongo(elements):
    entries = []
    unique_entries={}
    elements_with_H_O = elements
    elements_with_H_O.append('H')
    elements_with_H_O.append('O')
    elements_with_H_O.sort()
    print(elements)
    allcombinations = getAllCombinations(elements_with_H_O)
    myclient = pymongo.MongoClient("mongodb://localhost:27017/")
    mydb = myclient["mp"]
    mycol = mydb["aml_all5"]
    for c in allcombinations:
        c.sort()
        p='^'
        for i in c:
            p=p+i+'[0-9]+\s*' 
        p=p+'$'
        mq = { "elements" : c }
        md = mycol.find(mq)
        for x in md:
            print(x['pretty_formula'])
            fe=x['formation_energy_per_atom']
            if fe is not None:
                try:
                    fe=float(x['formation_energy_per_atom'])
                    if unique_entries.get(x['pretty_formula']) is not None:
                        if unique_entries[x['pretty_formula']] > x['formation_energy_per_atom']*x['natoms']:
                            unique_entries[x['pretty_formula']]=x['formation_energy_per_atom']*x['natoms']
                    else:
                        unique_entries[x['pretty_formula']]=x['formation_energy_per_atom']*x['natoms']
                except ValueError:
                    fe=100.0
                    print(x['task_id'], x['pretty_formula'], fe)
    for u in unique_entries:
        print(u, unique_entries[u])
        ent = PDEntry(Composition(u), unique_entries[u], attribute={'task_id':x['task_id']})
        ent = PourbaixEntry(ent)
        #ent.phase_type='Solid'
        entries.append(ent)
    return entries

def pourbaix(comp,target=None, range_ph=[0,16], range_voltage=[-2,2]):
    entries = []    
    try: 
        elements = [str(e) for e in np.unique(comp.split('-')).tolist()]
        if 'O' in elements: elements.remove('O')
        if 'H' in elements: elements.remove('H')
        if len(elements) == 5: print('2.5 milions etries will be considered. It will take more than 2 hours.')
        '''
        elif len(elements) > 6: 
            print('Too many etries needed to be considered. Please contact APOLEX admin (yongwoo.s@samsung.com) to perform this calculation.')
            return time.strftime("%c"), -2, [0,20], [-3,3]
        '''
    except ValueError as e:
        #print('Wrong input for the chemical space. Please put it with format "A-B-C-D"')
        #return time.strftime("%c"), -1, [0,20], [-3,3]
        return time.strftime(str(e)), -1, [0,20], [-3,3]
    elements.sort()
    elem = ''
    for i in range(len(elements)):
        elem = elem+"{}".format(elements[i])
        if i+1 == len(elements): elem = elem+''
        else: elem = elem+'-'
    print(elem)

    # Need to modify the directory
    printtimestamp('1')
    #fname = "/home/kimsooil/venv3/pourbaix/data/data_cache2/{}".format(elem)
    fname = "data/data_cache2/{}".format(elem)
    if os.path.isfile(fname) == True:
        with open(fname, 'rb') as fr: entries = pickle.load(fr)
        print('Found it-1')
    else:
        try:
            mpr = MPRester("8kEI0dFbzgfxX5cT")
            #entries = mpr.get_pourbaix_entries(elements)
            entries = getPourbaixEntryfromMongo(elements)
            #print(entries)
            with open(fname,'wb') as fw: pickle.dump(entries, fw, protocol=pickle.HIGHEST_PROTOCOL)
        except ValueError as e:
            if len(elements)>4:
                #return time.strftime("%c"), -2,  [0,20], [-3,3]
                return time.strftime(str(e)), -2,  [0,20], [-3,3]
            else:
                #return time.strftime("%c"), -1,  [0,20], [-3,3]
                return time.strftime(str(e)), -1,  [0,20], [-3,3]

    # Need to modify the directory
    #pbxname = "/home/kimsooil/venv3/pourbaix/data/data_cache2/{}-pbx".format(elem)
    pbxname = "data/data_cache2/{}-pbx".format(elem)
        
    printtimestamp('2')
    
    if os.path.isfile(pbxname) == True: 
        with open(pbxname, 'rb') as pbxr: pbx = pickle.load(pbxr)
        print('Found it-2')
    else:
        printtimestamp('2.5')
        pbx = PourbaixDiagram(entries,filter_solids=True)
        with open(pbxname,'wb') as pbxw: pickle.dump(pbx, pbxw, protocol=pickle.HIGHEST_PROTOCOL)
        
    printtimestamp('3')
    
    plotter = PourbaixPlotter(pbx)
    plt1 = plotter.get_pourbaix_plot(limits=[range_ph, range_voltage])
    
    printtimestamp('4')
    
    plt1.ylabel("E (V) vs SHE")
    #plt1.title("Pourbaix diagram of {} system".format(elem+'-O-H'),  fontsize=35, y=1.08)
    #plt1.tight_layout()
    #time_stamp = time.strftime("%c")
    time_stamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S.%f")
    plt1.savefig('static/pourbaix1_' + time_stamp + '.png')
    #print(entries)
    if target != None:    
        try:
            #entry = [entry for entry in entries if Composition(entry.composition).reduced_formula == Composition(target).reduced_formula][0]#942733-LLZO Li7La3Zr2O12
            entry = [entry for entry in entries if entry.name == target][0]
            print(entry)
            plt2 = plotter.plot_entry_stability(entry, pH_range=range_ph, pH_resolution=100,V_range=range_voltage, V_resolution=300, e_hull_max=1,limits=[range_ph, range_voltage])
            ##plt2.tight_layout()
            #plt2.title("Stability of {} in Pourbaix diagram of {} system".format(target, elem),  fontsize=30, y=1.08)
            plt2.ylabel("E (V) vs SHE")
            plt2.savefig('static/pourbaix2_' + time_stamp + '.png')
            images = 2
        except:
            images = 1
    else: images = 1
    return time_stamp, images, range_ph, range_voltage

#time_stamp, no_images, range_ph, range_voltage = pourbaix(comp='Li-Hf-La', target='La2Hf2O7', range_ph=[0,16], range_voltage=[-2,2])

