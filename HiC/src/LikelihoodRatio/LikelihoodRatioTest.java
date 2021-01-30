package LikelihoodRatio;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import mixedmodel.CompoundAnalyzer;
import mixedmodel.MultiPhenotype;

public class LikelihoodRatioTest {
	
	
	public static void permutatoin(String input_pheno_file, String output_pheno_file, int rounds_of_permutation, int index) {
		try {
			BufferedReader input_phenoReader = new BufferedReader(new FileReader(input_pheno_file));
			BufferedWriter output_phenoWriter = new BufferedWriter(new FileWriter(output_pheno_file));
			//read target pheno into arraylist
			ArrayList<String> tag_pheno_list =  new ArrayList<String>();
			ArrayList<String> tag_pheno_ID =  new ArrayList<String>();
			int tag_pheno_index=index+1; //as the 0 indicates ID, 1 indicates the 1st pheno
			String tag_pheno_name = input_phenoReader.readLine().split("\t")[tag_pheno_index];
			String tag_pheno_line = input_phenoReader.readLine();
			while(tag_pheno_line!=null) {
				String[] tag_pheno_lines = tag_pheno_line.split("\t");
				tag_pheno_ID.add(tag_pheno_lines[0]);
				tag_pheno_list.add(tag_pheno_lines[tag_pheno_index]);
				tag_pheno_line = input_phenoReader.readLine();
			}
			//generate permuated phenotype [col_0 for the real pheno, col_1 for the 1st permutated phenot, etc.]
			String output_pheno_header="#ID\t"+tag_pheno_name+"_real";
			for(int rp=0;rp<rounds_of_permutation;rp++) { //ID pheno_real pheno_rp1 pheno_rp2
				output_pheno_header=output_pheno_header+"\t"+tag_pheno_name+"_rp"+Integer.toString(rp+1);
			}
			output_phenoWriter.write(output_pheno_header+"\n");
			
			int sample_size = tag_pheno_ID.size();
			Random rand = new Random();
			for(int i=0;i<sample_size;i++) {
				String output_pheno_line=tag_pheno_ID.get(i)+"\t"+tag_pheno_list.get(i); //GTEx_1 0.45(pheno_real)
				for(int rp=0;rp<rounds_of_permutation;rp++) {
					int random_index = rand.nextInt(sample_size); //generate a random index from 0 to sample_size-1
					while(random_index==i) {random_index = rand.nextInt(sample_size);}
					output_pheno_line=output_pheno_line+"\t"+tag_pheno_list.get(random_index);
				}
				output_phenoWriter.write(output_pheno_line+"\n");
			}
			input_phenoReader.close();
			output_phenoWriter.close();
		}catch (Exception e) {e.printStackTrace();}
	}
	
	public static void LLRT(String genotype_hdf5, String input_pheno, String out_phe_file, String compound_file, int index, String global_kinship_file) {
		CompoundAnalyzer local_k=new CompoundAnalyzer(genotype_hdf5, compound_file);
		MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
		boolean plot=false;
		int min_sample_size=40;
		local_k.compound_VO(phenotypeS.phenotypes[index], genotype_hdf5, global_kinship_file, out_phe_file, plot, min_sample_size);
	}
	
	public static void LLRT_p(String phenotype_folder_br, String phenotype_folder_sr, String pheno_name, String compound_name,String output_pvalue_file, int rounds_of_permutation) {

		File f_br = new File(phenotype_folder_br);
		File f_sr = new File(phenotype_folder_sr);
		// This filter will only include files ending with .r.csv
//		FilenameFilter filter = new FilenameFilter() {
//		        @Override
//		        public boolean accept(File f, String name) {
//		            return name.endsWith(".r.csv");
//		        }
//		};
//		pathnames_br=f_br.list(filter);
//		pathnames_sr=f_sr.list(filter);

		ArrayList<String> pathnames_br_list = new ArrayList<>();
		ArrayList<String> pathnames_sr_list = new ArrayList<>();
		File[] f_br_listOfFiles = f_br.listFiles();
		File[] f_sr_listOfFiles = f_sr.listFiles();
		for (int i = 0; i < f_br_listOfFiles.length; i++) {
			  if (f_br_listOfFiles[i].isFile()) {
				  pathnames_br_list.add(f_br_listOfFiles[i].getName());
			  }
		}
		for (int i = 0; i < f_sr_listOfFiles.length; i++) {
			  if (f_sr_listOfFiles[i].isFile()) {
				  pathnames_sr_list.add(f_sr_listOfFiles[i].getName());
			  }
		}
		System.out.println(pathnames_br_list.get(0));
		System.out.println(pathnames_br_list.size());
		System.out.println(pathnames_sr_list.size());
		
		if(pathnames_br_list.size()==rounds_of_permutation+1 && pathnames_sr_list.size()==rounds_of_permutation+1) {
			System.out.println("Finishing analysis for all rounds of permtuations");
		}else {
			String real="Compound_VO.0."+pheno_name+"_real.r.csv";
			if(!pathnames_br_list.contains(real)) {System.out.println("Didn't finish analyzing phenotype round br real");}
			if(!pathnames_sr_list.contains(real)) {System.out.println("Didn't finish analyzing phenotype round sr real");}
			for(int i=1;i<rounds_of_permutation+1;i++) {
				//Compound_VO.320.ENSG00000149527.13_PLCH2_rp320.r.csv
				String tag="Compound_VO."+i+"."+pheno_name+"_rp"+i+".r.csv";
				if(!pathnames_br_list.contains(tag)) {System.out.println("Didn't finish analyzing phenotype round br "+i);}
				if(!pathnames_sr_list.contains(tag)) {System.out.println("Didn't finish analyzing phenotype round sr "+i);}
			}
			System.exit(0);
		}
		
		//read compounds information
		 try {
			 String[] compound_name_arr = compound_name.split("_");
			 String compound_r1 = compound_name_arr[0]+"_1:1:1";
			 String compound_r2 = compound_name_arr[1]+"_1:1:1";
			 BufferedWriter one_compound_raw_bWriter =  new BufferedWriter(new FileWriter(output_pvalue_file+".raw"));
			 BufferedWriter one_compound_LL_bWriter =  new BufferedWriter(new FileWriter(output_pvalue_file+".LL_LLR"));
			 one_compound_raw_bWriter.write("#round,compound,best_ml,p-value,variance_local,variance_global,variance_e\n");
			 one_compound_LL_bWriter.write("#round,LL_int,LL_r1,LL_r2,LLR_1,LLR_2\n");
			 double LLR_1_pvalue=-1,LLR_2_pvalue=-1; //real
			 ArrayList<Double> LLR_p_1_list= new ArrayList<Double>();
			 ArrayList<Double> LLR_p_2_list= new ArrayList<Double>();
			 //get real LL_int,LL_r1,LL_r2,LLR_1,LLR_2
			 String cp_br_file=phenotype_folder_br+"/Compound_VO.0."+pheno_name+"_real.r.csv";
			 String cp_sr_file=phenotype_folder_sr+"/Compound_VO.0."+pheno_name+"_real.r.csv";
			 String resString = Get_LL_LLR("real", cp_br_file, cp_sr_file, compound_name, compound_r1, compound_r2, one_compound_LL_bWriter, one_compound_raw_bWriter);
			 String[] resStrings = resString.split(",");
			 double LLR_1 = Double.parseDouble(resStrings[4]);
			 double LLR_2 = Double.parseDouble(resStrings[5]);
			 //permutated LL_int,LL_r1,LL_r2,LLR_1,LLR_2
			 for(int rp=1;rp<rounds_of_permutation+1;rp++) {
				 String cp_p_br_file=phenotype_folder_br+"/Compound_VO."+rp+"."+pheno_name+"_rp"+rp+".r.csv";
				 String cp_p_sr_file=phenotype_folder_sr+"/Compound_VO."+rp+"."+pheno_name+"_rp"+rp+".r.csv";
				 String res_p_String = Get_LL_LLR("r"+rp, cp_p_br_file, cp_p_sr_file, compound_name, compound_r1, compound_r2, one_compound_LL_bWriter, one_compound_raw_bWriter);
				 String[] res_p_Strings = res_p_String.split(",");
				 double LLR_p_1 = Double.parseDouble(res_p_Strings[4]);
				 double LLR_p_2 = Double.parseDouble(res_p_Strings[5]);
				 if(!Double.isInfinite(LLR_p_1)&&!Double.isNaN(LLR_p_1)) LLR_p_1_list.add(LLR_p_1);
				 if(!Double.isInfinite(LLR_p_2)&&!Double.isNaN(LLR_p_2)) LLR_p_2_list.add(LLR_p_2);
			 }
			 one_compound_raw_bWriter.close();
			 one_compound_LL_bWriter.close();
			 Collections.sort(LLR_p_1_list);
			 Collections.sort(LLR_p_2_list);
			 System.out.println("Number of permutation for r1 "+LLR_p_1_list.size());
			 System.out.println("Number of permutation for r2 "+LLR_p_2_list.size());
			 for(int i=0;i<LLR_p_1_list.size();i++) {
				 if(Double.compare(LLR_1, LLR_p_1_list.get(i))>0)continue;
				 else {LLR_1_pvalue=(i+1.0)/LLR_p_1_list.size();}
			 }
			 for(int i=0;i<LLR_p_2_list.size();i++) {
				 if(Double.compare(LLR_2, LLR_p_2_list.get(i))>0)continue;
				 else {LLR_2_pvalue=(i+1.0)/LLR_p_2_list.size();}
			 }
			 System.out.println(pheno_name+","+compound_name+",pvalue_r1,"+LLR_1_pvalue+",pvalue_r2,"+LLR_2_pvalue+"\n");			 
		 }catch (Exception e) {e.printStackTrace();}
	}
	
	public static String Get_LL_LLR(String index, String cp_br_file, String cp_sr_file, String compound_name,String compound_r1, String compound_r2, BufferedWriter LLWriter, BufferedWriter rawWriter) {
		String resString=""; 
		try {
			 double LL_int=0.0,LL_r1=0.0,LL_r2=0.0,LLR_1=0.0,LLR_2=0.0;	 
			 BufferedReader cp_brReader = new BufferedReader(new FileReader(cp_br_file));
			 BufferedReader cp_srReader = new BufferedReader(new FileReader(cp_sr_file));
			 String cp_br_line = cp_brReader.readLine();
			 String cp_sr_line = cp_srReader.readLine();
			 while(cp_br_line!=null) {
				 String[] cp_br_line_arr = cp_br_line.split(",");
				 if(cp_br_line_arr[0].equals(compound_name)) {
					 LL_int=Double.parseDouble(cp_br_line_arr[1]);
					 rawWriter.write(index+","+cp_br_line+"\n");
					 break;
				 }
				 cp_br_line=cp_brReader.readLine();
			 }cp_brReader.close();
			 boolean r1_written=false,r2_written=false;
			 while(cp_sr_line!=null&&(!r1_written||!r2_written)) {
				 String[] cp_sr_line_arr = cp_sr_line.split(",");
				 if(cp_sr_line_arr[0].equals(compound_r1)) {
					 LL_r1=Double.parseDouble(cp_sr_line_arr[1]);
					 rawWriter.write(index+","+cp_sr_line+"\n");
					 r1_written=true;
				 }else if(cp_sr_line_arr[0].equals(compound_r2)) {
					 LL_r2=Double.parseDouble(cp_sr_line_arr[1]);
					 rawWriter.write(index+","+cp_sr_line+"\n");
					 r1_written=true;
				 }
				 cp_sr_line=cp_srReader.readLine();
			 }
			 cp_srReader.close();
			 if(r1_written) {LLR_1=LL_int-LL_r1;}
			 else {LL_r1=Double.NaN;LLR_1=Double.NaN;}
			 if(r2_written) {LLR_2=LL_int-LL_r2;}
			 else {LL_r2=Double.NaN;LLR_2=Double.NaN;}
			 resString=index+","+LL_int+","+LL_r1+","+LL_r2+","+LLR_1+","+LLR_2;
			 LLWriter.write(resString+"\n");
		}catch (Exception e) {e.printStackTrace();}
		return resString;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			System.out.println("======================================================\n"
					+ "Developer: Qing Li Dec 2020\n"
					+ "Usage: java -Xmx4g -jar LLR.jar [function] <compulsory parameters> [optional parameters]\n"
					+ "Supported functions:\n"
					+ "\tpermutation\n"
					+ "\tLLRT (LikelihoodRatioTest)\n"
					+ "\tLLRT_p (pvalue for likelihoodratio)");
			System.exit(0);
		}else{
			String function = args[0];
			if(function.equals("permutation")) {
				int number_of_parameter=3;
				if(args.length<number_of_parameter*2+1) {
					System.out.println("Usage: \n"
							+ "\t<-ipf\tphenotype>\n"
							+ "\t<-p\trounds of permutation>\n"
							+ "\t<-opf\toutput phenotype file name>\n"
							+ "\t[-index\tindex of pheotype to be permutated (df=0 for pheno file contains 2 columns[ID,Pheno])]\n"
							+ "There are "+number_of_parameter+" mandatory paraneters. \nBut you have only specified "
							+(args.length-1)/2+" parameters");
				}
				else {
					String input_pheno_file=null;
					String output_pheno_file=null;
					int rounds_of_permutation=0;
					int index=0;
					for(int k=1;k<args.length;k++) {
						if(args[k].startsWith("-")) {
							if(args[k].startsWith("-ipf"))input_pheno_file=args[k+1];
							else if(args[k].startsWith("-p"))rounds_of_permutation=Integer.parseInt(args[k+1]);
							else if(args[k].startsWith("-opf"))output_pheno_file=args[k+1];
							else if(args[k].startsWith("-index"))index=Integer.parseInt(args[k+1]);
						}
					}
					if(input_pheno_file==null||output_pheno_file==null||rounds_of_permutation==0) {
						System.out.println("input or output can't be null!");
					}else {
						permutatoin(input_pheno_file, output_pheno_file, rounds_of_permutation,index);
					}
				}
					
			}else if(function.equals("LLRT")) {
				int number_of_parameter=5;
				if(args.length<number_of_parameter*2+1) {
					System.out.println("Usage: \n"
							+ "\t<-ig\tgenotype (hdf5)>\n"
							+ "\t<-ip\tphenotype>\n"
							+ "\t<-ik_g\tglobal kinship>\n"
							+ "\t<-of\toutpout folder> [end with /]\n"
							+ "\t<-ic\tcompound_file>\n"
							+ "\t[-index\tindex of phenotype requires analysis] (0 for phenotype file has 2 cols[IDs,Pheno])\n"
							+ "There are "+number_of_parameter+" mandatory paraneters. \nBut you have only specified "
							+(args.length-1)/2+" parameters");
				}
				else {
					String genotype_hdf5=null;
					String input_pheno=null;
					String output_folder=null;
					String compound_file=null;
					String global_kinship=null;
					int the_phe_index=-1;
					for(int k=1;k<args.length;k++){
						if(args[k].startsWith("-")){
							if(args[k].equals("-ig"))genotype_hdf5=args[k+1];
							else if(args[k].equals("-ip"))input_pheno=args[k+1];
							else if(args[k].equals("-ik_g"))global_kinship=args[k+1];
							else if(args[k].equals("-of"))output_folder=args[k+1];
							else if(args[k].equals("-ic"))compound_file=args[k+1];
							else if(args[k].equals("-index"))the_phe_index=Integer.parseInt(args[k+1]);
						}
					}if(genotype_hdf5==null||input_pheno==null||global_kinship==null||output_folder==null||compound_file==null) {
						System.out.println("input or output can't be null!");
					}else {
						MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
						if(the_phe_index==-1){
							System.out.println("Running all phenotypes? It is suggested to specify a phenotype index. " +
									"Otherwise it may be slow." +
									"\nLet us have a try!");
							for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++) {
								String out_phe_file=output_folder+"Compound_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".r.csv";
								LLRT(genotype_hdf5, input_pheno, out_phe_file, compound_file, phe_index, global_kinship);}
							}else {
								String out_phe_file=output_folder+"Compound_VO."+the_phe_index+"."+phenotypeS.phenotypes[the_phe_index].phe_id+".r.csv";
								LLRT(genotype_hdf5, input_pheno, out_phe_file, compound_file, the_phe_index, global_kinship);
							}
					}
				}
				
			}else if(function.equals("LLRT_p")) {
				int number_of_parameter=6;
				if(args.length<number_of_parameter*2+1) {
					System.out.println("Usage: \n"
							+ "\t<-io_br\tinput folder (all results for phenoytpe with double regions)>\n"
							+ "\t<-io_sr\tinput folder (all results for phenoytpe with single region)>\n"
							+ "\t<-ic\tone compound name for double regions>\n"
							+ "\t<-pheno\tphenotype name>\n"
							+ "\t<-p\trounds of permutation>\n"
							+ "\t<-op_file\toutput pvalue file name>\n"
							+ "There are "+number_of_parameter+" mandatory paraneters. \nBut you have only specified "
							+(args.length-1)/2+" parameters");
				}
				else {
					String phenotype_folder_br=null;
					String phenotype_folder_sr=null;
					String compound_name=null;
					String output_pvalue_file=null;
					String pheno_name=null;
					int rounds_of_permutation=0;
					
					for(int k=1;k<args.length;k++){
						if(args[k].startsWith("-")){
							if(args[k].equals("-io_br"))phenotype_folder_br=args[k+1];
							else if(args[k].equals("-io_sr"))phenotype_folder_sr=args[k+1];
							else if(args[k].equals("-ic"))compound_name=args[k+1];
							else if(args[k].equals("-pheno"))pheno_name=args[k+1];
							else if(args[k].equals("-p"))rounds_of_permutation=Integer.parseInt(args[k+1]);
							else if(args[k].equals("-op_file"))output_pvalue_file=args[k+1];
						}
					}if(phenotype_folder_br==null||phenotype_folder_sr==null||compound_name==null||output_pvalue_file==null||pheno_name==null) {
						System.out.println("input or output can't be null!");
					}else {
						LLRT_p(phenotype_folder_br, phenotype_folder_sr, pheno_name, compound_name, output_pvalue_file, rounds_of_permutation);
					}
					}
			}
		
		}
			
	}



}
