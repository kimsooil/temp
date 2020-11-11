@app.route('/organic_compound_search_cep/', methods=['GET'])
#@app.cache.cached(timeout=0)
def organic_compound_search_cep():
    result={}
    output_data= []

    smiles = request.args.get('smiles')
    inchi = request.args.get('inchi') # no inchi is considered for CEP organic db
    formula = request.args.get('formula')
    id = request.args.get('id')
    min_gap = request.args.get('min_gap')
    max_gap = request.args.get('max_gap')
    print(smiles, formula,inchi,id)

    ranges = request.args.get('ranges')
    '''
    if (smiles is "" or smiles is None) and (inchi is "" or inchi is None) and (formula is "" or formula is None) and (id is "" or id is None) and (min_gap is "" or min_gap is None) and (max_gap is "" or max_gap is None):
        return jsonify({"SMILES":smiles, "InChI":inchi, "formula":formula, "id":id, "output":[]})
    '''
    mq={}
    if smiles is not None and smiles is not "":
        #the following replacements are needed only when re.escape() is not used
        #smiles=smiles.replace('(', '\(').replace(')', '\)').replace('[', '\[').replace(']', '\]').replace('#', '\#').replace('=', '\=').replace('+', '\+') 
        #mq["SMILES_str"]={"$regex":smiles+"\s"}
        mq["SMILES_str"]=re.compile('^' + re.escape(smiles) + '$', re.IGNORECASE)
    '''
    if inchi is not None and inchi is not "":
        inchi=inchi.replace('(', '\(').replace(')', '\)').replace('[', '\[').replace(']', '\]').replace('#', '\#').replace('=', '\=').replace('+', '\+')
        mq["InChI"]={"$regex":inchi+"\s"}
    '''
    if formula is not None and formula is not "":
        mq['stoich_str']={"$regex":formula}
    if id is not None and id is not "":
        id=int(id)
        #mq['id']={"$regex":id}
        mq['id']=id

    if ranges is not None and ranges is not "":
        ranges = [float(x) for x in ranges.split('_')]
        mq['e_gap_alpha']={'$gte':ranges[0], '$lte':ranges[1]}
        mq['e_homo_alpha']={'$gte':ranges[2], '$lte':ranges[3]}
        mq['e_lumo_alpha']={'$gte':ranges[4], '$lte':ranges[5]}

    print(mq)
    #md= (myclient["organic"])["xyz_gdb9"].find(mq)
    md= (myclient["organic"])["cep"].find(mq).limit(1000) #get first 1000 results only....
    for x in md:
        out_temp = json.dumps(x,default=str)
        out_temp = json.loads(out_temp)
        output_data.append(out_temp)

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S.%f")

    result={'timestamp':timestamp, 'SMILES':smiles, 'formula':formula, 'output':output_data, 'length':len(output_data)}

    return jsonify(result)
    
@app.route('/organic_compound_entry_cep/', methods=['GET'])
#@app.cache.cached(timeout=0)
def organic_compound_entry_cep():
    result={}
    output_data= []

    id = request.args.get('id')

    if id is "" or id is None:
        return jsonify({"id":id, "output":[]})

    mq={}
    if id is not None and id is not "":
        #id=re.sub(r'[^a-zA-Z0-9_\.]', '', id)
        id=re.sub(r'[^0-9_\.]', '', id)
        id=int(id)
        mq={"id":id}

    print(mq)
    #md= (myclient["organic"])["xyz_gdb9"].find(mq)
    md= (myclient["organic"])["cep"].find(mq)

    for x in md:
        out_temp = json.dumps(x,default=str)
        out_temp = json.loads(out_temp)
        output_data.append(out_temp)

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S.%f")

    result={'timestamp':timestamp, 'id':id, 'output':output_data, 'length':len(output_data)}

    return jsonify(result)
