@app.route('/icsd_anonsearch/<string:anonformula>/', methods=['GET'])
#@app.cache.cached(timeout=0)
def icsd_anonsearch(anonformula='ABCD4'):
    include = request.args.get('include') # eg. /anonsearch/ABCD4/?include=Li,Fe,P,O/
    exclude = request.args.get('exclude')
    sg_sym = request.args.get('sg_sym')
    sg_num = request.args.get('sg_num')
    crystal_system=request.args.get('crystal_system')

    elements=[]
    output_entries=[]

    alpha_pattern=[]
    if exclude:
        ex_elements=exclude.split(',')
        for ex_ele in ex_elements:
            alpha_pattern.append({"alphabetical_formula":{"$not":re.compile(ex_ele+"[0-9]+")}})
    if include is not None:
        # check ntypes first (how many different elements)
        ntypes_in_anonformula = sum(c.isalpha() for c in anonformula)
        #print('ntypes_in_anonformula', ntypes_in_anonformula)
        natoms_in_anonformula = 0 # if * used in anonformula, natoms is not considered. Only ntypes is considered
        if anonformula.find("*")==-1: # not contain '*'
            numbers = re.findall(r'[a-zA-Z]([\d+[\.\d+]?]?)', anonformula)
            #print('numbers', numbers)
            for number in numbers:
                if number=='':
                    natoms_in_anonformula += 1 
                else:
                    natoms_in_anonformula += int(number)
            #print("natoms in anonformula", natoms_in_anonformula)
            elements=include.split(',')
            elements.sort()
            for element in elements:
                alpha_pattern.append({"alphabetical_formula":{"$regex":element+"[0-9]+"}})
            if ntypes_in_anonformula > len(elements):
                mq= {"anonymized_formula":anonformula, "ntypes":ntypes_in_anonformula, "$and":alpha_pattern}
            else:
                mq= {"anonymized_formula":anonformula, "$and":alpha_pattern}

            if sg_sym is not None and sg_sym!='' and sg_sym!='null':
                alpha_pattern.append({"spacegroup_short":sg_sym})
            if sg_num is not None and sg_num!='' and sg_num!='null':
                alpha_pattern.append({"sym_num":int(sg_num)})
            if crystal_system is not None and crystal_system!='' and crystal_system!='null':
                crystal_system = crystal_system.capitalize()
                alpha_pattern.append({"crystal_system":crystal_system})
        else: # wildcard(*) in anonymized_formula
            elements=include.split(',')
            elements.sort()
            for element in elements:
                alpha_pattern.append({"alphabetical_formula":{"$regex":element+"[0-9]+"}})
            mq= {"ntypes":ntypes_in_anonformula, "$and":alpha_pattern}

        md= (myclient["icsd"])["icsd_data4"].find(mq)

        for x in md:
            out_temp = json.dumps(x, default=str)
            out_temp = json.loads(out_temp)
            output_entries.append(out_temp)

    else: # just generic (eg. ABC2, ABCD4) without 'include'
        #print("alpha_pattern", alpha_pattern)
        ntypes_in_anonformula = sum(c.isalpha() for c in anonformula)
        mq ={"anonymized_formula":str(anonformula), "$and":alpha_pattern}

        md =  (myclient["icsd"])["icsd_data4"].find(mq)

        for x in md:
            out_temp = json.dumps(x, default=str)
            out_temp = json.loads(out_temp)
            output_entries.append(out_temp)

    result={'anonymized_formula':anonformula, 'include_elements':include, 'exclude_elements':exclude, 'spacegroup_symbol':sg_sym, 'spacegroup_number':sg_num, 'crystal_system':crystal_system, 'output':output_entries, 'length':len(output_entries)}

    return jsonify(result)
