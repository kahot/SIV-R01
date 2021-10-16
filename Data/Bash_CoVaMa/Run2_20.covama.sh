
#1. run ViReMa.py 

mkdir Run2_20
python2.7 ~/programs/ViReMa_0.7_4CoVaMa/ViReMa.py SIV_R01/SIVB670 ~/programs/SIV-R01/Output/PID_Consensus/Run2_20.consensus.fasta Run2_20/Run2_20 -F --Output_Tag Run2_20 --Output_Dir Run2_20/ --p 8


#2 Run CoVaMa_Make_Matrices.py 
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Make_Matrices.py Run2_20/Run2_20 SIV/AY032751env_padded.fasta  NT --ViReMa_Output1 Run2_20_CoVaMa_Output.txt --Min_Fusion_Coverage 10 --Mode2 Nucs
  
  
#3. Run CoVaMa_Analyse_Matrices.py
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Analyse_Matrices.py Run2_20/Run2_20.Total_Matrices.py.pi Output/Run2_20.Output.txt Run2_20/ --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted