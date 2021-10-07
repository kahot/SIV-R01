
#1. run ViReMa.py 

mkdir Run2_9_
python2.7 ~/programs/ViReMa_0.7_4CoVaMa/ViReMa.py SIV_R01/SIVB670 ~/programs/SIV-R01/Output/PID_Consensus/Run2_9_.consensus.fasta Run2_9_/Run2_9_ -F --Output_Tag Run2_9_ --Output_Dir Run2_9_/ --p 8


#2 Run CoVaMa_Make_Matrices.py 
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Make_Matrices.py Run2_9_/Run2_9_ SIV/AY032751env_padded.fasta  NT --ViReMa_Output1 Run2_9__CoVaMa_Output.txt --Min_Fusion_Coverage 10 --Mode2 Nucs
  
  
#3. Run CoVaMa_Analyse_Matrices.py
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Analyse_Matrices.py Run2_9_/Run2_9_.Total_Matrices.py.pi Output/Run2_9_.Output.txt Run2_9_/ --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted
