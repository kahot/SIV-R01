
#1. run ViReMa.py 

mkdir Run5_12
python2.7 ~/programs/ViReMa_0.7_4CoVaMa/ViReMa.py SIV_R01/SIVB670 ~/programs/SIV-R01/Output/PID_Consensus/Run5_12.consensus.fasta Run5_12/Run5_12 -F --Output_Tag Run5_12 --Output_Dir Run5_12/ --p 8


#2 Run CoVaMa_Make_Matrices.py 
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Make_Matrices.py Run5_12/Run5_12 SIV/AY032751env_padded.fasta  NT --ViReMa_Output1 Run5_12_CoVaMa_Output.txt --Min_Fusion_Coverage 10 --Mode2 Nucs
  
  
#3. Run CoVaMa_Analyse_Matrices.py
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Analyse_Matrices.py Run5_12/Run5_12.Total_Matrices.py.pi Output/Run5_12.Output.txt Run5_12/ --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted
