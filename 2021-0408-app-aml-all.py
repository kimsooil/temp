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
        #composition=re.sub(r'\W+', '', composition)
        comp=None
        try:
            comp = Composition(composition)
        except:
            pass
        elements = [str(e) for e in comp.elements]
        elements.sort()
        #elements=result['elements']
        for el in ['Li', 'Na', 'Mg', 'K', 'Ca', 'Al', 'Zn']:
            if el in elements and 'O' in elements:
                evolution=get_voltage_profile(composition)
    result['evolution']=evolution

    mpid=taskid.replace('aml-m-', 'mp-')
    activation_energy_barrier=(myclient["softBV"])["activation_energy_barrier"].find({"id":mpid})
    e_a=[]
    if activation_energy_barrier is not None:
        for x in activation_energy_barrier:
            a=json.dumps(x, default=str)
            e_a.append(json.loads(a))
    result['activation_energy_barrier']=e_a

    return jsonify(result)
  
  @app.route('/icsd_element/<string:composition>/', methods=['GET'])
#@app.cache.cached(timeout=0)
def icsd_element(composition='Li-La-Zr-O'):
    result={}
    entries = []
    unique_entries = {} 
    outdata=[]
    composition_entered=composition
    composition=re.sub(r'\W+', '', composition_entered)
    comp=None
    try:
        comp=Composition(composition)
    except:
        return jsonify({"composition":composition, "output":[]})
    elements = [str(e) for e in comp.elements]
    elements.sort()
    allcombinations = getAllCombinations(elements)

    mydb_pd = myclient["icsd"]
    mycol_pd = mydb_pd["icsd_data4"]

    for c in allcombinations:
        c.sort()
        p='^'
        for i in c:
            p=p+i+'[0-9]+\s*' 
        p=p+'$'
        mq = { "elements" : c } # elements should be sorted alphabetically
        md = mycol_pd.find(mq)
        for x in md:
            out_temp = json.dumps(x,default=str)
            out_temp = json.loads(out_temp)
            outdata.append(out_temp)
            
            #fe=get_formation_energy_from_aml(x['icsd_id']) # SLOWWWWWW
            fe=0.0
            unique_entries[x['composition']]=fe
            print(x['icsd_id'], x['composition'], fe)
            
            if unique_entries.get(x['composition']) is not None:
                if unique_entries[x['composition']] > fe:
                    unique_entries[x['composition']]=fe
            
    for u in unique_entries:
        print(u, unique_entries[u])
        #entries.append(PDEntry(u, unique_entries[u], attribute={'task_id':x['task_id']}))
        entries.append(PDEntry(Composition(u), unique_entries[u], attribute={'icsd_id':x['icsd_id']}))
    
    print(entries)
    timestamp=""
    pd_filename=""
    lines=[]
    lines2=[]
    stable_entries=[]
    unstable_entries=[]
    stable=[]
    unstable=[]
    plotter=None
    if len(elements)<5:
        pd = PhaseDiagram(entries)

        try:
            #plotter=PDPlotter(pd)
            plotter=PDPlotter(pd, show_unstable=1)
        except:
            return jsonify({'composition':composition, 'output':[], 'pd_file_name':''})
        lines, stable_entries,unstable_entries=plotter.pd_plot_data
        #print(lines)
        #print(stable_entries)

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S.%f")
        pd_filename="pd-"+timestamp+".svg"
        plotter.write_image("static/"+pd_filename)

        stable=[]
        for k, v in stable_entries.items():
            #print(k, v)
            s=[]
            p=str(k).replace('(','').replace(')','').split(', ')
            s.append(float(p[0]))
            s.append(float(p[1]))
            if len(elements)==4: s.append(float(p[2]))
            s.append(str(v).replace('PDEntry :', '').replace(' with energy = ', ':'))
            stable.append(s)
        unstable=[]
        for k, v in unstable_entries.items():
            #print('key=',k,'value=',v)
            u=[]
            p=str(v).replace('(','').replace(')','').split(', ')
            u.append(float(p[0]))
            u.append(float(p[1]))
            if len(elements)==4: u.append(float(p[2]))
            u.append(str(k).replace('PDEntry :', '').replace(' with energy = ', ':'))
            unstable.append(u)
        '''
        lines2=[]
        for i in range(len(lines)):
            aline=[]
            aline1=[]
            aline2=[]
            #print(i, lines[i][0][0], lines[i][0][1], lines[i][1][0], lines[i][1][1])
            aline1.append(lines[i][0][0]) #x1
            aline1.append(lines[i][1][0]) #y1
            if len(elements)==4: aline1.append(lines[i][2][0])
            aline2.append(lines[i][0][1]) #x2
            aline2.append(lines[i][1][1]) #y2
            if len(elements)==4: aline2.append(lines[i][2][1])
            aline.append(aline1)
            aline.append(aline2)
            lines2.append(aline)
        '''
    lines = json.dumps(lines, cls=NumpyArrayEncoder)
    
    url = "http://0.0.0.0:7000/pd2/"+composition+"/" # endining with / is must. it's related to flask's url format
    response = urlopen(url)
    mp_pd_data = json.loads(response.read())

    lines2=mp_pd_data['lines'] # use pd lines of mp element search

    pd_plot_data=""
    if plotter!=None:
        pd_plot_data=str(plotter.pd_plot_data)
    
    result={'composition':composition_entered, 'elements':elements, 'timestamp':timestamp, 'pd_filename':pd_filename, "unique_entries":unique_entries, 'lines':lines2, 'stable':stable, 'unstable':unstable, 'pd_plot_data':pd_plot_data,'output':outdata, 'length':len(outdata), "mp_pd":mp_pd_data, "stable_length":len(stable), "unstable_length":len(unstable), "unique_length":len(unique_entries)}
   # result={'composition':composition_entered, 'elements':elements, 'output':outdata, 'length':len(outdata)}

    return jsonify(result)
