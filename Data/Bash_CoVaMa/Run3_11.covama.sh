
#1. run ViReMa.py 

mkdir Run3_11
python2.7 ~/programs/ViReMa_0.7_4CoVaMa/ViReMa.py SIV_R01/SIVB670 ~/programs/SIV-R01/Output/PID_Consensus/Run3_11.consensus.fasta Run3_11/Run3_11 -F --Output_Tag Run3_11 --Output_Dir Run3_11/ --p 8


#2 Run CoVaMa_Make_Matrices.py 
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Make_Matrices.py Run3_11/Run3_11 SIV/AY032751env_padded.fasta  NT --ViReMa_Output1 Run3_11_CoVaMa_Output.txt --Min_Fusion_Coverage 10 --Mode2 Nucs
  
  
#3. Run CoVaMa_Analyse_Matrices.py
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Analyse_Matrices.py Run3_11/Run3_11.Total_Matrices.py.pi Output/Run3_11.Output.txt Run3_11/ --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted