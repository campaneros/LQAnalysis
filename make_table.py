import argparse

def parse_file(filename):
    data = {}
    with open(filename, 'r') as file:
        for line in file:
            key, value = line.split(':')
            key = key.split('"')[1]  # extract the string within quotes
            data[key] = float(value.split(" ")[0])  # convert the value to float
    return data

def main(file1, file2, file3):
    # Parse both files
    data1 = parse_file(file1)
    data2 = parse_file(file2)
    data3 = parse_file(file3)

    # Create a new dictionary to hold the results
    result = {}
    result2 = {}

    # Divide the corresponding values
    for key in data1:
        if key in data2:
            result[key] = data2[key] / data1[key]
    for key in data2:
        for key in data3:
            result2[key] = data3[key] / data2[key]

    # Print the results
    for key, value in result.items():
        print('"' + key + '":' + str(value))

    # Write the result  to a file
    filename=file1.replace('.txt', '') 
    app=filename.split('_')
    print(app)
    index=app.index('limits')
    print(app[index])
    if "bmu" in app:
        filename = app[index]+"_"+app[index+1]+"_"+app[index+2]+"_result.txt"
    else:  
        filename = app[index]+"_"+app[index+1]+"_"+app[index+2]+"_result.txt"
    with open(filename, 'w') as file:
        for key, value in result.items():
            file.write('4cat/2cat L '+str(app[index+2])+" " + key + '":' + str(value) + " 2cat: "+str(data1[key])+" 4cat "+str(data2[key])+'\n')
        for key, value in result2.items():
            file.write('8cat/4cat L '+str(app[index+2])+" " + key + '":' + str(value) + " 4cat: "+ str(data2[key])+" 8cat "+str(data3[key])+'\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process three input files.")
    parser.add_argument('file1', help='First input file')
    parser.add_argument('file2', help='Second input file')
    parser.add_argument('file3', help='Second input file')
    args = parser.parse_args()
    main(args.file1, args.file2,args.file3)
