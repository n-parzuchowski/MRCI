names= ["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al", \
        "Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe", \
        "Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr", \
        "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn"] 

a = raw_input("Enter Nucleus (C14): ")

for I in range(len(a)):
    if a[I].isdigit():
        elem = a[:I]
        ATarg = int(a[I:])
        ZTarg = names.index(elem)
        NTarg = ATarg- ZTarg
        break

a = raw_input("Enter Reference (C14): ")


for I in range(len(a)):
    if a[I].isdigit():
        elem = a[:I]
        Aref = int(a[I:])
        Zref = names.index(elem)
        Nref = Aref- Zref
        break

yes = raw_input("is this a single reference state? (y/n): " )

yes = yes.lower()

par = raw_input("Enter desired parity (0,1): ")

inter  = raw_input("Enter interaction tag: (chi2b3b_srg0625): ")

a = raw_input("Enter comma delimited eMax: ")
eMaxes = a.strip().split(",")

a = raw_input("Enter comma delimited hw: ")
hws = a.strip().split(",")

reffile = raw_input("Enter reference file name: ")


a = raw_input("enter comma delimited MRCI walltimes for each eMax (just an integer for hours): ")
timesMRCI = a.strip().split(",")

a = raw_input("enter comma delimited MRCI memories for each eMax (just an integer for gb): ")
memsMRCI = a.strip().split(",")

a = raw_input("enter comma delimited IMSRG memories for each eMax (just an integer for gb): ")
memsIMSRG = a.strip().split(",")

fq = open("submit_mrci.bat","w")
fq.write("#!/bin/bash \n\n")


fx = open("submit_imsrg.bat","w")
fx.write("#!/bin/bash \n\n")

fxx = open("submit_restart.bat","w")
fxx.write("#!/bin/bash \n\n")

fy = open("submit_hfb.bat","w")
fy.write("#!/bin/bash \n\n")


fz = open("submit_nord.bat","w")
fz.write("#!/bin/bash \n\n")

etick =0
for e in eMaxes:
    for hw in hws:
    ### get prefix
        if len(e) == 1:
            ex = "0"+e
        else:
            ex = e 
        if len(hw) == 2:
            hwx = "0"+hw
        else:
            hwx = hw

        if (int(e)> 10):
            prefix = elem+str(ATarg)+"_"+inter+"_eMax"+ex+"lMax10_hwHO"+hwx
            prefix2 = elem+str(Aref)+"_"+inter+"_eMax"+ex+"lMax10_hwHO"+hwx
        else:
            prefix = elem+str(ATarg)+"_"+inter+"_eMax"+ex+"_hwHO"+hwx
            prefix2 = elem+str(Aref)+"_"+inter+"_eMax"+ex+"_hwHO"+hwx

        intfile = prefix2 + ".ham.me2b.gz" 
        rhofile = prefix2 + ".lambda.me2b.gz"
        spfile = "hk"+e+".sps"


        ### write MRCI inifile

        fl = open(prefix+".ini","w")

        fl.write("### Enter number of nucleons (Z,N) in reference nucleus\n")
        fl.write(str(Zref) + "  "+str(Nref)+"\n") 
        fl.write("### Enter number of nucleons (Z,N) in target nucleus\n")
        fl.write(str(ZTarg) + "  "+str(NTarg)+"\n")        
        fl.write("### Enter parity (0,1)\n")
        fl.write(par+"\n") 
        fl.write("### Enter sp file\n")
        fl.write(spfile+"\n")
        fl.write("### Enter ham.me2b file\n")
        fl.write(intfile+"\n")
        fl.write("### Enter lambda.me2b file\n")
        fl.write(rhofile+"\n")
        fl.write("### Enter .ref file\n")
        fl.write(reffile+"\n")

        fl.close()
        
        ### write MRCI pbs file

        time = timesMRCI[etick] 
        time = "0"+time+":00:00" 
        mem = memsMRCI[etick]+"gb"
        ppn = "16"
        
        fl =open("pbsMRCI_"+prefix,"w")
        
        fl.write("#!/bin/sh\n\n")
        fl.write("#PBS -l walltime="+time+"\n")
        fl.write("#PBS -l nodes=1:ppn="+ppn+"\n")
        fl.write("#PBS -l mem="+mem+"\n")
        fl.write("#PBS -j oe\n")
        fl.write("#PBS -N "+prefix+"_MRCI\n")
        fl.write("#PBS -M nathan.parz@gmail.com\n")
        fl.write("#PBS -m a\n\n")
        fl.write("cd $HOME/MRCI/src\n\n")
        fl.write("export OMP_NUM_THREADS="+ppn+"\n\n")
        fl.write("./run_mrci "+prefix+".ini\n\n")
        fl.write("qstat -f ${PBS_JOBID}\nexit 0")

        fl.close()

        fq.write("qsub pbsMRCI_"+prefix+"\n")


        ### write HFB pbs file

        time = "04:00:00"
        ppn = "16"    
        mem = "10gb"
        
        fl=open("pbsHFB_"+prefix2,"w")
        
        fl.write("#!/bin/sh\n\n")
        fl.write("#PBS -l walltime="+time+"\n")
        fl.write("#PBS -l nodes=1:ppn="+ppn+"\n")
        fl.write("#PBS -l mem="+mem+"\n")
        fl.write("#PBS -j oe\n")
        fl.write("#PBS -N "+prefix2+"_HFB\n")
        fl.write("#PBS -M nathan.parz@gmail.com\n")
        fl.write("#PBS -m a\n\n")
        fl.write("cd $HOME/hfb\n")
        fl.write("source to_source_bfr_make\n\n")
        fl.write("export OMP_NUM_THREADS="+ppn+"\n\n")
        fl.write("./solve_hfb Nucl="+names[Zref]+str(Aref)+" eMax="+e \
                 +" E3Max=14 l3Max=10 hwHO="+hw+" IntID="+inter+"\n\n")
        fl.write("qstat -f ${PBS_JOBID}\nexit 0")

        fl.close()

        fy.write("qsub pbsHFB_"+prefix2+"\n")


        if (int(e) > 10 ):
            thing = "_lMax10"
        else:
            thing = ""
            
        ### write NORD pbs file

        
        time = "04:00:00"
        ppn = "16"    
        mem = "10gb"        
        
        fl=open("pbsNORD_"+prefix2,"w")
        
        fl.write("#!/bin/sh\n\n")
        fl.write("#PBS -l walltime="+time+"\n")
        fl.write("#PBS -l nodes=1:ppn="+ppn+"\n")
        fl.write("#PBS -l mem="+mem+"\n")
        fl.write("#PBS -j oe\n")
        fl.write("#PBS -N "+prefix2+"_NORD\n")
        fl.write("#PBS -M nathan.parz@gmail.com\n")
        fl.write("#PBS -m a\n\n")
        fl.write("cd $HOME/nsuite\n")
        fl.write("source to_source_bfr_make\n\n")
        fl.write("export OMP_NUM_THREADS="+ppn+"\n\n")
        if (yes == "n"):
            fl.write("./normalorder_ham RefType=PNP res/"+names[Zref]+str(Aref)+\
                     "_eMax"+ex+thing+"_hwHO"+hwx+".hfbc\n\n")
        else:
            fl.write("./normalorder_ham res/"+names[Zref]+str(Aref)+"_eMax"+ex \
                     +thing+"_hwHO"+hwx+".hfbc\n\n")

        fl.write("qstat -f ${PBS_JOBID}\nexit 0")

        fl.close()

        fz.write("qsub pbsNORD_"+prefix2+"\n")



        ### write FLINT pbs file

        
        time = "24:00:00"
        mem = memsIMSRG[etick]+"gb"
        ppn = "16"
        
        
        fl=open("pbsFLINT_"+prefix2,"w")
        
        fl.write("#!/bin/sh\n\n")
        fl.write("#PBS -l walltime="+time+"\n")
        fl.write("#PBS -l nodes=1:ppn="+ppn+"\n")
        fl.write("#PBS -l mem="+mem+"\n")
        fl.write("#PBS -j oe\n")
        fl.write("#PBS -N "+prefix2+"_FLINT\n")
        fl.write("#PBS -M nathan.parz@gmail.com\n")
        fl.write("#PBS -m a\n\n")
        fl.write("cd $HOME/nsuite\n")
        fl.write("source to_source_bfr_make\n\n")
        fl.write("export OMP_NUM_THREADS="+ppn+"\n\n")

        if (yes == "n"):
            fl.write("./solveimsrg_flint RefType=PNP EtaType=BrillouinMinimal sMax=5000.0 res/" \
                     +names[Zref]+str(Aref)+"_eMax"+ex+thing+"_hwHO"+hwx+".hfbc\n\n")
        else:
            fl.write("./solveimsrg_flint EtaType=BrillouinMinimal sMax=5000.0 res/" \
                     +names[Zref]+str(Aref)+"_eMax"+ex+thing+"_hwHO"+hwx+".hfbc\n\n") 

        fl.write("qstat -f ${PBS_JOBID}\nexit 0")

        fl.close()

        fx.write("qsub pbsFLINT_"+prefix2+"\n")


                ### write FLINT pbs file

        
        time = "24:00:00"
        mem = memsIMSRG[etick]+"gb"
        ppn = "16"
        
        
        fl=open("pbsRESTART_"+prefix2,"w")
        
        fl.write("#!/bin/sh\n\n")
        fl.write("#PBS -l walltime="+time+"\n")
        fl.write("#PBS -l nodes=1:ppn="+ppn+"\n")
        fl.write("#PBS -l mem="+mem+"\n")
        fl.write("#PBS -j oe\n")
        fl.write("#PBS -N "+prefix2+"_RESTART\n")
        fl.write("#PBS -M nathan.parz@gmail.com\n")
        fl.write("#PBS -m a\n\n")
        fl.write("cd $HOME/nsuite\n")
        fl.write("source to_source_bfr_make\n\n")
        fl.write("export OMP_NUM_THREADS="+ppn+"\n\n")

        fl.write("cd res \n")
        fl.write("ln -s ${SCRATCH}/"+inter+"/res/"+prefix2+".flint.chk \n\n" )
        fl.write("./solveimsrg_flint RefType=PNP EtaType=BrillouinMinimal sMax=5000.0 res/" \
                 +prefix2+".flint.chk\n\n")

        fl.write("qstat -f ${PBS_JOBID}\nexit 0")

        fl.close()

        fxx.write("qsub pbsRESTART_"+prefix2+"\n")

    etick+=1
fx.close()
fy.close()
fz.close()
fq.close()
