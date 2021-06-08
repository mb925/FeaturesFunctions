def parse_obo(obo_file):
    rel_list = []
    go_terms = {} #dictionary key = go term id value dict containing info on the go term (name index....)
    # temp variables used to remember what go term we are parsing 
    temp_id = ''
    temp_namespace = '' 
    temp_name = ''
    temp_def = ''
    temp_alt_id = []
    with open(obo_file) as f:
        for line in f:
            line = line.strip().split(": ")
            if line and len(line) > 1:
                k, v = line[:2]
                if k == "id" and (v.startswith("GO:") or v.startswith("HP:") or v.startswith("DO:")):
                    # when a new id is found first we have to input the entry in the go term dict
                    # also we will have to input a new entry when we reach EOF 
                    if temp_id != '':
                        go_terms[temp_id] = {'name': temp_name, 'namespace': temp_namespace, 'def': temp_def, 'parents': [], 'children': [], 'alt_ids': temp_alt_id}
                        temp_alt_id = []
                    temp_id = v
                elif k == "alt_id" and (v.startswith("GO:") or v.startswith("HP:") or v.startswith("DO:")):
                    temp_alt_id.append(v)
                elif k == "name":
                    temp_name = v
                elif k == "namespace" and v != 'external':
                    temp_namespace = v
                elif k == "def":
                    temp_def = v
                elif k == "is_a" and (v.startswith("GO:") or v.startswith("HP:") or v.startswith("DO:")):
                    # add (temp_id,v) tuple to the relation list
                    s = v.split('!')[0].strip()
                    rel_list.append([temp_id, s, temp_namespace])
                """elif k == "relationship" and v.startswith("part_of"):
                    s = v.split()[1].strip()
                    rel_list.append([temp_id, s, temp_namespace])"""
              

        go_terms[temp_id] = {'name': temp_name, 'namespace': temp_namespace, 'def': temp_def, 'parents': [], 'children': [], 'alt_ids': temp_alt_id}             

    for r in rel_list:
        i = r[0]
        j = r[1]
        go_terms[i]['parents'].append(j)
        go_terms[j]['children'].append(i)
    return go_terms