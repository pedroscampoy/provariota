list_res_position = []
dict_res_position_base = {}

with open("MTB_Base_Calibration_List.vcf", 'r') as f:
    for line in f:
        if line.startswith("#"):
            next(f)
        else:
            list_line = line.split("\t")
            position = list_line[1]
            nucleotide = list_line[4].upper()
            resistance = list_line[-1].split(";RES=")[-1]
            list_res_position.append(position)
            dict_res_position_base[position] = [nucleotide, resistance] 
            #print(position, nucleotide, resistance)

with open("MTB_Resistance_Mediating.txt", 'r') as f:
    for _ in range(1):
        next(f)
    for line in f:
        list_line = line.split("\t")
        if list_line[2] == "SNP":
            position = list_line[1]
            nucleotide = list_line[5].upper()
            resistance = list_line[21]
            if position not in list_res_position and not resistance.startswith("phylo") :
                list_res_position.append(position)
                dict_res_position_base[position] = [nucleotide, resistance]

print(dict_res_position_base)