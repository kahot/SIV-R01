
#1. run ViReMa.py 

mkdir Run5_6_
python2.7 ~/programs/ViReMa_0.7_4CoVaMa/ViReMa.py SIV_R01/SIVB670 ~/programs/SIV-R01/Output/PID_Consensus/Run5_6_.consensus.fasta Run5_6_/Run5_6_ -F --Output_Tag Run5_6_ --Output_Dir Run5_6_/ --p 8


#2 Run CoVaMa_Make_Matrices.py 
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Make_Matrices.py Run5_6_/Run5_6_ SIV/AY032751env_padded.fasta  NT --ViReMa_Output1 Run5_6__CoVaMa_Output.txt --Min_Fusion_Coverage 10 --Mode2 Nucs
  
  
#3. Run CoVaMa_Analyse_Matrices.py
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Analyse_Matrices.py Run5_6_/Run5_6_.Total_Matrices.py.pi Output/Run5_6_.Output.txt Run5_6_/ --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted
