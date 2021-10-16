
#1. run ViReMa.py 

mkdir Run0_8_
python2.7 ~/programs/ViReMa_0.7_4CoVaMa/ViReMa.py SIV_R01/SIVB670 ~/programs/SIV-R01/Output/PID_Consensus/Run0_8_.consensus.fasta Run0_8_/Run0_8_ -F --Output_Tag Run0_8_ --Output_Dir Run0_8_/ --p 8


#2 Run CoVaMa_Make_Matrices.py 
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Make_Matrices.py Run0_8_/Run0_8_ SIV/AY032751env_padded.fasta  NT --ViReMa_Output1 Run0_8__CoVaMa_Output.txt --Min_Fusion_Coverage 10 --Mode2 Nucs
  
  
#3. Run CoVaMa_Analyse_Matrices.py
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Analyse_Matrices.py Run0_8_/Run0_8_.Total_Matrices.py.pi Output/Run0_8_.Output.txt Run0_8_/ --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted