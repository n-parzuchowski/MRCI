import sys
import os.path 

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


print "Target Z,N: ",ZTarg,NTarg
a = raw_input("Enter Reference (C14): ")


for I in range(len(a)):
    if a[I].isdigit():
        elem = a[:I]
        Aref = int(a[I:])
        Zref = names.index(elem)
        Nref = Aref- Zref
        break

print "Reference Z,N: ",Zref,Nref
yes = raw_input("is this a single reference state? (y/n): " )

yes = yes.lower()

par = raw_input("Enter desired parity (0,1): ")

inter  = raw_input("Enter interaction tag: (chi2b3b_srg0625): ")
if (inter ==""):
    inter= "chi2b3b_srg0625"

SCR= "/mnt/ls15/scratch/users/parzuch6"
hfb_resdir = "/mnt/home/parzuch6/hfb/res/"+inter
scratch_resdir = SCR+"/"+inter+"/res"

short_resdir = "res/"+inter

nsuite_resdir = "/mnt/home/parzuch6/nsuite/res/"+inter

if not os.path.isdir(hfb_resdir):
    os.system("mkdir "+hfb_resdir)

if not os.path.isdir(nsuite_resdir):
    os.system("mkdir "+nsuite_resdir)

if not os.path.isdir(scratch_resdir):
    if not os.path.isdir(SCR+"/"+inter):
        os.system("mkdir "+SCR+"/"+inter )
    os.system("mkdir "+scratch_resdir)

a = raw_input("Enter comma delimited eMax: ")
eMaxes = a.strip().split(",")

a = raw_input("Enter comma delimited hw: ")
hws = a.strip().split(",")

reffile = raw_input("Enter reference file name: ")


a = raw_input("enter comma delimited MRCI walltimes for each eMax (just an integer for hours): ")
timesMRCI = a.strip().split(",")

a = raw_input("enter comma delimited MRCI memories for each eMax (just an integer for gb): ")
memsMRCI = a.strip().split(",")

a = raw_input("should we make IMSRG input files too? (y/n): ")
if (a.lower() == "y"):
    imsrg = True
else:
    imsrg = False

if imsrg:
    a = raw_input("enter comma delimited IMSRG memories for each eMax (just an integer for gb): ")
    memsIMSRG = a.strip().split(",")


fq = open(names[ZTarg]+str(ATarg)+"_"+inter+"_mrci.bat","w")
fq.write("#!/bin/bash \n\n")

if imsrg:
    fx = open(names[Zref]+str(Aref)+"_"+inter+"_imsrg.bat","w")
    fx.write("#!/bin/bash \n\n")

    fxx = open(names[Zref]+str(Aref)+"_"+inter+"_restart.bat","w")
    fxx.write("#!/bin/bash \n\n")
    
    fy = open(names[Zref]+str(Aref)+"_"+inter+"_hfb.bat","w")
    fy.write("#!/bin/bash \n\n")
    
    fz = open(names[Zref]+str(Aref)+"_"+inter+"_nord.bat","w")
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
            prefix = names[ZTarg]+str(ATarg)+"_"+inter+"_eMax"+ex+"_lMax10_hwHO"+hwx
            prefix2 = names[Zref]+str(Aref)+"_"+inter+"_eMax"+ex+"_lMax10_hwHO"+hwx
        else:
            prefix = names[ZTarg]+str(ATarg)+"_"+inter+"_eMax"+ex+"_hwHO"+hwx
            prefix2 = names[Zref]+str(Aref)+"_"+inter+"_eMax"+ex+"_hwHO"+hwx

        intfile = prefix2 + ".ham.me2b.gz" 
        rhofile = prefix2 + ".lambda.me2b.gz"

        if (e < 12):
            spfile = "hk"+e+".sps"
        else:
            spfile = "hk"+e+"_lmax10.sps"
        
        ### write MRCI inifile

        fl = open("inifiles/"+prefix+".ini","w")

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


        if not imsrg:
            continue
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

        if (int(e) > 10):
            fl.write("./solve_hfb Nucl="+names[Zref]+str(Aref)+" eMax="+e \
                         +" E3Max=14 lMax=10 hwHO="+hw+" IntID="+inter+\
                         " ResDir="+hfb_resdir+"\n\n")
        else:
            fl.write("./solve_hfb Nucl="+names[Zref]+str(Aref)+" eMax="+e \
                         +" E3Max=14 hwHO="+hw+" IntID="+inter+\
                         " ResDir="+hfb_resdir+"\n\n")
            
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

        mem = "14gb"        
        
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
        fl.write("cp "+hfb_resdir+"/"+names[Zref]+str(Aref)+\
                     "_eMax"+ex+thing+"_hwHO"+hwx+".hfbc "+nsuite_resdir+"\n\n")
        if (yes == "n"):
            fl.write("./normalorder_ham RefType=PNP ResDir="+short_resdir+\
                         " "+short_resdir+"/"+names[Zref]+str(Aref)+\
                         "_eMax"+ex+thing+"_hwHO"+hwx+".hfbc\n")
            fl.write("./normalorder_obs RefType=PNP ResDir="+short_resdir+\
                         " IntID=Hcm_beta01.00 "+short_resdir+"/"+names[Zref]+str(Aref)+\
                         "_eMax"+ex+thing+"_hwHO"+hwx+".hfbc\n\n")
        else:
            fl.write("./normalorder_ham  ResDir="+short_resdir+\
                         " "+short_resdir+"/"+names[Zref]+str(Aref)+"_eMax"+ex \
                         +thing+"_hwHO"+hwx+".hfbc\n")
            fl.write("./normalorder_obs  ResDir="+nsuite_resdir+\
                         " IntID=Hcm_bet01.00 "+short_resdir+"/"+names[Zref]+str(Aref)+"_eMax"+ex \
                         +thing+"_hwHO"+hwx+".hfbc\n\n")

        fl.write("rm "+nsuite_resdir+"/"+names[Zref]+str(Aref)+\
                     "_eMax"+ex+thing+"_hwHO"+hwx+".hfbc\n\n")
        fl.write("qstat -f ${PBS_JOBID}\nexit 0")

        fl.close()

        fz.write("qsub pbsNORD_"+prefix2+"\n")



        ### write MAGNUS pbs file

        
        time = "24:00:00"
        mem = memsIMSRG[etick]+"gb"
        ppn = "16"
        
        
        fl=open("pbsMAGNUS_"+prefix2,"w")
        
        fl.write("#!/bin/sh\n\n")
        fl.write("#PBS -l walltime="+time+"\n")
        fl.write("#PBS -l nodes=1:ppn="+ppn+"\n")
        fl.write("#PBS -l mem="+mem+"\n")
        fl.write("#PBS -j oe\n")
        fl.write("#PBS -N "+prefix2+"_MAGNUS\n")
        fl.write("#PBS -M nathan.parz@gmail.com\n")
        fl.write("#PBS -m a\n\n")
        fl.write("cd $HOME/nsuite\n")
        fl.write("source to_source_bfr_make\n\n")
        fl.write("export OMP_NUM_THREADS="+ppn+"\n\n")

        if (yes == "n"):
            fl.write("./solveimsrg_magnus ObsID=Hcm_beta01.00 RefType=Multi ResDir="\
                         +short_resdir+" EtaType=BrillouinMinimal WriteAll sMax=5000.0 "+short_resdir \
                         +"/"+prefix2+".mref\n\n")
        else:
            fl.write("./solveimsrg_magnus ObsID=Hcm_beta01.00 ResDir="\
                         +short_resdir+" EtaType=BrillouinMinimal WriteAll sMax=5000.0 "+short_resdir \
                         +"/"+prefix2+".mref\n\n")
            
        fl.write("qstat -f ${PBS_JOBID}\nexit 0")

        fl.close()

        fx.write("qsub pbsMAGNUS_"+prefix2+"\n")

        
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

        fl.write("cd  "+nsuite_resdir+"\n")
        fl.write("ln -s ${SCRATCH}/"+inter+"/res/"+prefix2+".magnus.chk \n" )
        fl.write("cd $HOME/nsuite \n\n")
        fl.write("./solveimsrg_magnus RefType=Multi ResDir="\
                     +short_resdir+" EtaType=BrillouinMinimal WriteAll sMax=5000.0 " \
                     +short_resdir+"/"+prefix2+".magnus.chk\n\n")

        fl.write("qstat -f ${PBS_JOBID}\nexit 0")

        fl.close()

        fxx.write("qsub pbsRESTART_"+prefix2+"\n")

    etick+=1


fq.close()

os.system("chmod 0755 "+names[ZTarg]+str(ATarg)+"_"+inter+"_mrci.bat")

if imsrg:
    fx.close()
    fxx.close()
    fy.close()
    fz.close()
    os.system("chmod 0755 "+names[Zref]+str(Aref)+"_"+inter+"_imsrg.bat")
    os.system("chmod 0755 "+names[Zref]+str(Aref)+"_"+inter+"_restart.bat")
    os.system("chmod 0755 "+names[Zref]+str(Aref)+"_"+inter+"_hfb.bat")
    os.system("chmod 0755 "+names[Zref]+str(Aref)+"_"+inter+"_nord.bat")
