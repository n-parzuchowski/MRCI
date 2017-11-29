import glob
import os,sys

def get_id_num(fname):
    ind = fname.index(".o")
    return int(fname[ind+2:])

files = []
for fname in glob.iglob("*.o*"):
    # check if fname is actually an output file
    
    try: 
        id_num = get_id_num(fname)
        if (id_num > 50000):
          files.append(fname)
        else: 
            print(fname)
            sys.exit()
    except ValueError:
        print("Excluding file: "+fname)    
        
files=sorted(files)

for f in files:
    if ("MAGNUS" in f): 
        fl=open(f,"r") 
        f_contents=fl.readlines()
        fl.close()
        
        clean = True 
        for row in f_contents:
            if "killed" in row.lower():
                print("JOB KILLED")
                print(f) 
                print(row)
                ### okay we need to figure out if this is the most recent file or not 
                cur_id = get_id_num(f)
                nextfile = f[:f.index("MAGNUS")]+"RESTART"
                for g in files: 
                    if nextfile in g:
                        ### found same system 
                        ### check if it is a more a recent run 
                        next_id = get_id_num(g)      
                        if (next_id > cur_id):
                            os.system("rm "+f) 
                            files.remove(f)
                            break
                        else:
                            # file is older
                            print("\n Removing older file:")
                            print(g) 
                            os.system("rm "+g) 
                            files.remove(g)
                clean = False
                break
            if "fault" in row.lower(): 
                if "fault_tolerant" in row.lower(): continue                   
                print("SEG FAULT") 
                print(f)
                print(row) 
                clean = False
                break 

        if clean:
            print("file "+f+" indicates convergence in MAGNUS run") 


for f in files:
    if ("RESTART" in f): 
        fl=open(f,"r") 
        f_contents=fl.readlines()
        fl.close()
        
        clean = True 
        for row in f_contents:
            if "killed" in row.lower():
                print("JOB KILLED")
                print(f) 
                print(row)
                ### okay we need to figure out if this is the most recent file or not 
                cur_id = get_id_num(f)
                nextfile = f[:f.index(".o4")]
                for g in files: 
                    if (nextfile in g) and g != f:
                        ### found same system 
                        ### check if it is a more a recent run 
                        next_id = get_id_num(g)      
                        if (next_id > cur_id):
                            os.system("rm "+f) 
                            files.remove(f)
                            break
                        elif (next_id < cur_id):
                            # file is older
                            print("\n Removing older file:")
                            print(g) 
                            os.system("rm "+g) 
                            files.remove(g)
                        else: 
                            print("something is wrong")

                clean = False
                break
            if "fault" in row.lower(): 
                if "fault_tolerant" in row.lower(): continue                   
                print("SEG FAULT") 
                print(f)
                print(row) 
                clean = False
                break 

        if clean:
            print("file "+f+" indicates convergence in RESTART run") 
        
            
        
