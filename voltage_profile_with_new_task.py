@app.route('/task/<string:taskid>/', methods=['GET'])
#@app.cache.cached(timeout=0)
def task_info(taskid='aml-m-1222978'):
    result={}
    task_found =  (myclient["mp"])["aml_all5"].find_one({"task_id":taskid})
    if task_found is not None:
        result= json.dumps(task_found, default=str)
        result= json.loads(result)
    else:
        task_found=(myclient["mp"])["aml_all5"].find_one({"task_ids":{"$regex":taskid}}) # exact match
        if task_found is not None:
            result= json.dumps(task_found, default=str)
            result= json.loads(result)
    evolution=[]
    if result['pretty_formula'] is not None:
        composition=result['pretty_formula']
        elements=result['elements']
        for el in ['Li', 'Na', 'Mg', 'K', 'Ca', 'Al', 'Zn']:
            if el in elements and 'O' in elements:
                evolution=get_voltage_profile(composition)
    result['evolution']=evolution
    return jsonify(result)
    
def get_voltage_profile(composition="LiCoO2", chargeElement="Li"):
    unique_entries={}
    entries = []
    composition=re.sub(r'\W+', '', composition)
    comp=None
    try:
        comp = Composition(composition)
    except:
        return jsonify({"composition":composition, "output":[]})
    elements = [str(e) for e in comp.elements]
    elements.sort()

    allcombinations = getAllCombinations(elements)

    mydb_pd = myclient["mp"]
    mycol_pd = mydb_pd["aml_all5"]

    for c in allcombinations:
        c.sort()
        p='^'
        for i in c:
            p=p+i+'[0-9]+\s*' 
        p=p+'$'
        mq = { "elements" : c }
        md = mycol_pd.find(mq)
        for x in md:
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
    for u in unique_entries:
        entries.append(PDEntry(u, unique_entries[u], attribute={'task_id':x['task_id']}))

    evolution=[]
    if len(elements)<5:
        try:
            pd=PhaseDiagram(entries)
            #plotter=PDPlotter(pd, show_unstable=False)
            evol=pd.get_element_profile(chargeElement,comp)
            element_energy = evol[0]["chempot"]

            for i, d in enumerate(evol):
                evolution.append(str(-(d["chempot"] - element_energy))+":"+str(d['evolution'])+":"+str(d['reaction']))
                print(-(d["chempot"] - element_energy), d['evolution'], d['reaction'])
        except:
            return []
    return evolution

@app.route('/voltage_profile/<string:composition>/', methods=['GET'])
#@app.cache.cached(timeout=0)
def voltage_profile(composition='Li7La3Zr2O12', chargeElement="Li"):
    result={}
    evolution=get_voltage_profile(composition, chargeElement)
    result={'evol':evolution}
    return jsonify(result)

@app.route('/voltage_profile_mpr/<string:composition>/', methods=['GET'])
#@app.cache.cached(timeout=0)
def voltage_profile_mpr(composition='Li7La3Zr2O12'):
    result={}
    composition=re.sub(r'\W+', '', composition)
    comp=None
    try:
        comp = Composition(composition)
    except:
        return jsonify({"composition":composition, "output":[]})
    elements = [str(e) for e in comp.elements]
    elements.sort()

    evolution=[]
    if len(elements)<5:

        try:
            mpr=MPRester(MAPI_KEY)
            compat=MaterialsProjectCompatibility()
            unprocessed_entries=mpr.get_entries_in_chemsys(elements)
            processed_entries=compat.process_entries(unprocessed_entries)

            pd=PhaseDiagram(processed_entries)

            #plotter=PDPlotter(pd, show_unstable=False)
            evol=pd.get_element_profile("Li",comp)
            element_energy = evol[0]["chempot"]

            for i, d in enumerate(evol):
                evolution.append(str(-(d["chempot"] - element_energy))+":"+str(d['evolution'])+":"+str(d['reaction']))
                print(-(d["chempot"] - element_energy), d['evolution'], d['reaction'])
        except:
            return jsonify({'composition':composition, 'evol':[]})

    result={'evol':evolution}
    return jsonify(result)
