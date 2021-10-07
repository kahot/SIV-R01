
#1. run ViReMa.py 

mkdir Run6_5_
python2.7 ~/programs/ViReMa_0.7_4CoVaMa/ViReMa.py SIV_R01/SIVB670 ~/programs/SIV-R01/Output/PID_Consensus/Run6_5_.consensus.fasta Run6_5_/Run6_5_ -F --Output_Tag Run6_5_ --Output_Dir Run6_5_/ --p 8


#2 Run CoVaMa_Make_Matrices.py 
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Make_Matrices.py Run6_5_/Run6_5_ SIV/AY032751env_padded.fasta  NT --ViReMa_Output1 Run6_5__CoVaMa_Output.txt --Min_Fusion_Coverage 10 --Mode2 Nucs
  
  
#3. Run CoVaMa_Analyse_Matrices.py
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Analyse_Matrices.py Run6_5_/Run6_5_.Total_Matrices.py.pi Output/Run6_5_.Output.txt Run6_5_/ --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted
