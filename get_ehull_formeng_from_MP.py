def get_e_hull_and_formation_energy_from_MP(vasp_data):
    E_hull=-1000.0
    E_form=-1000.0
    with MPRester(MAPI_KEY) as mpr:
        ent_list = list(set(vasp_data.atomic_symbols))
        entries = mpr.get_entries_in_chemsys(ent_list)
        up_entry = vasp_data.get_computed_entry(inc_structure=True)
        compatibility = MaterialsProjectCompatibility()
        up_entry = compatibility.process_entry(up_entry)
        entries_all = compatibility.process_entries(entries + [up_entry])
        pd = PhaseDiagram(entries_all)
        E_form = pd.get_form_energy_per_atom(up_entry)
        decomp, e_above_hull = pd.get_decomp_and_e_above_hull(up_entry)
        E_hull=e_above_hull
        print('E_form', E_form)
        '''
        my_entry = vasp_data.get_computed_entry(inc_structure=False)
        E_hull =mpr.get_stability([my_entry])[0]['e_above_hull']
        '''
        print('E_hull', E_hull)

    return E_hull, E_form

@app.route('/xml_to_cif/<string:filename>/', methods=['GET'])
#@app.cache.cached(timeout=0)
def get_cif_from_xml(filename='test.cif'):

    result={}
    cif_found=""
    if filename is not None:
        try:
            url="http://127.0.0.1/uploadedfiles/"+filename
            print(url)
            xml_found=urlopen(url).read().decode('utf-8')
            #print(xml_found)
            f = open("static/"+filename, "w")
            f.write(xml_found)
            f.close()
            path="static/"+filename
            #print(path)
            cif_found, out_read, structure=convert_xml_to_cif(path) # convert uploaded xml to cif
            cif_str=""
            if cif_found:
                with open(path+'.cif', 'r') as file:
                    cif_str = file.read()
            volume = structure.volume
            density = structure.density

            finder = SpacegroupAnalyzer(out_read.final_structure)
            s_group_num = finder.get_space_group_number()
            s_group_sbl = finder.get_space_group_symbol()
            crystal_sys = finder.get_crystal_system()

            n_atom = out_read.final_structure.num_sites
            band_gap, cbm, vbm, Eg_direct = out_read.eigenvalue_band_properties
            e_final = out_read.final_energy
            n_site = out_read.final_structure.num_sites

            E_hull=-1000.0
            E_form=-1000.0
            capacity=0.0
            try:
                E_hull, E_form=get_e_hull_and_formation_energy_from_MP(out_read)
            except:
                pass
            #E_hull = mpr.get_stability([entry4E_hull])[0]['e_above_hull']

            compounds_found_list=""
            formula = out_read.final_structure.formula
            comp=Composition(formula)
            anonymized_formula=comp.anonymized_formula
            #sorted_formula=sorted(formula.split())
            #alphabetical_formula=' '.join(sorted_formula)
            alphabetical_formula=comp.alphabetical_formula
            mq={"alphabetical_formula":alphabetical_formula}
            md= (myclient["mp"])["aml_all5"].find(mq)

            #structure1 = CifParser(path+".cif", occupancy_tolerance=100).get_structures()[0]
            for x in md:
                compounds_found_list += str( x['task_id'])+ "_"
                structure2 = CifParser.from_string(str( x['cif']), occupancy_tolerance=100).get_structures()[0]
                if structures_matched(structure, structure2)==True:
                    compounds_found_list += "same "
                else:
                    compounds_found_list += "different "
            
        except FileNotFoundError:
            return jsonify(result)

        if cif_found:
            result={"cif_file":(filename+'.cif'), "cif_str":cif_str, "spg_num":s_group_num, "spg_sym":s_group_sbl, "crystal_sys":crystal_sys, 
                    "formula":formula, "alphabetical_formula":alphabetical_formula, 
                    "compounds_list":compounds_found_list, "anonymized_formula":anonymized_formula,
                    "natoms":n_atom, "nsites":n_site, "band_gap":band_gap, "e_final":e_final, "volume":volume, "density":density, "capacity":capacity,
                    "E_hull":E_hull, "E_form":E_form}
    return jsonify(result)

@app.route('/inserting_new_to_aml_m', methods=['POST'])
def inserting_new_to_aml_m():
    if request.method=='POST':
        data = request.json

        mydb = myclient["mp"]
        mycol = mydb["aml_all5"]        
        new_material={}

        #print(data)
        ##cif=data['cif']
        cif_str=""
        cif_file=data['cif_file']
        cif_path="static/"+cif_file
        if cif_file:
            with open(cif_path, 'r') as file:
                cif_str = file.read()        
        now=datetime.datetime.now()
        time_stamp = now.strftime("%Y-%m%d-%H%M-%S") #2020-06-10-30
        time_stamp2= now.strftime("%y%m%d%H%M%S") # 20061030
        new_material['time_inserted']=time_stamp
        new_material['memo']="added by upload_compare app"
        new_material['task_id']="aml-m-"+time_stamp2
        new_material['pretty_formula']=data['formula'].replace(" ", "")
        new_material['alphabetical_formula']=data['alphabetical_formula']
        new_material['anonymized_formula']=data['anonymized_formula']
        new_material['ntypes']=len(data['alphabetical_formula'].split())
        new_material['cif']=cif_str
        new_material['spacegroup']={"number":data['spg_num'], "symbol":data['spg_sym'], "crystal_system":data['crystal_sys']}
        new_material['natoms']=int(data['natoms'])
        new_material['final_energy_per_atom']=float(data['e_final'])/float(data['natoms'])
        new_material['band_gap']=float(data['band_gap'])
        new_material['volume']=float(data['volume'])
        new_material['density']=float(data['density'])
        new_material['capacity']=float(data['capacity'])

        # don't know how to get it for now. but need initializing to avoid error later
        new_material['e_above_hull']=float(data['e_above_hull'])
        new_material['formation_energy_per_atom']=float(data['formation_energy_per_atom'])
        new_material['icsd_ids']=""
        new_material['decomposes_to']={}
        #new_material['formation_energy_per_atom']=-1000000.0
        #new_material['formation_energy_per_atom2']=-1000000.0

        id = mycol.insert_one(new_material)

        '''
        filename='static/new_material-'+cif_file+'.txt'
        file = open(filename,"w") 
        file.write(cif) 
        file.close()         
        '''
        
        response = jsonify({'result':'Success-POST', 'taskid':new_material['task_id']})
        response.headers.add('Access-Control-Allow-Origin', '*')
        return response
    else:
        response= jsonify({'result':'ok'})
        response.headers.add('Access-Control-Allow-Origin', '*')
        return response
