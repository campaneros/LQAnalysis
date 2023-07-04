import optparse
import ROOT




usage = "usage: python plotLimits.py -i /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/limit/Res1ToRes2ToGluGlu -l 40926"

parser = optparse.OptionParser(usage)

parser.add_option("-d", "--inputdir", dest="inputdir",
                  help="input directory with limits")



(opt, args) = parser.parse_args()

if not opt.inputdir:   
    parser.error('input directory not provided')





file_list = []
for path, subdirs, files in os.walk(opt.inputdir):
	for name in files:
		print(os.path.join(path, name))
		file_list.append(os.path.join(path, name))	

def parse_file(filename):
    data = {}
    with open(filename, 'r') as file:
        for line in file:
            key, value = line.split('||')
            key = key.split('"')[1]  # extract the string within quotes
            data[key] = float(value)  # convert the value to float
    return data