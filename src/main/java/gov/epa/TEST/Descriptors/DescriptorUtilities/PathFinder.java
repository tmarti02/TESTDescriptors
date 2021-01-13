package gov.epa.TEST.Descriptors.DescriptorUtilities;

import java.text.Collator;

import java.util.Collections;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;;

public class PathFinder {
    private static final Logger logger = LogManager.getLogger(PathFinder.class);
	
	
	public static void CalculateVertexDistanceDegrees(IAtomContainer m,int Distance[][],int [] vdd) {
				
		for (int i=0;i<=m.getAtomCount()-1;i++) {			
			vdd[i]=0;
			for (int j=0;j<=m.getAtomCount()-1;j++) {// use lower triangular distance submatrix
				vdd[i]+=Distance[i][j];
			}
		}
	}
	
	
	
	
	public static void CalculateDistanceMatrix(IAtomContainer m,int [][] A) {
		// this class was rewritten from computeFloydAPSP in
		// PathTools.java in CDK so that the distance matrix is not
		// created with new- so that can use same distance matrix in different
		// places
		int[][] C = org.openscience.cdk.graph.matrix.AdjacencyMatrix.getMatrix(m);
		
				
//		for (int f = 0; f < m.getBondCount(); f++){
//			IBond bond = m.getBond(f);
//	        int indexAtom1 = m.getAtomNumber(bond.getAtom(0));
//			int indexAtom2 = m.getAtomNumber(bond.getAtom(1));
//			System.out.println(f+"\t"+indexAtom1+"\t"+indexAtom2+"\t"+bond.getOrder());
//		}
		
		int i;
		int j;
		int k;
		int n = C.length;
		//int[][] A = new int[n][n];
		//System.out.println("Matrix size: " + n);
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (C[i][j] == 0) {
					A[i][j] = 999999999;
				}
				else {
					A[i][j] = 1;
				}
			}
		}
		for (i = 0; i < n; i++) {
			A[i][i] = 0;
			// no self cycle
		}
		for (k = 0; k < n; k++) {
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (A[i][k] + A[k][j] < A[i][j]) {
						A[i][j] = A[i][k] + A[k][j];
						//P[i][j] = k;        // k is included in the shortest path
					}
				}
			}
		}
		
		
	}
	
	
	
	
	public static LinkedList FindClusters(IAtomContainer m,int M) {

		//finds clusters of length 3 or 4 when have atom with 3 or 4 atoms attached
		
		
		LinkedList MasterList = new LinkedList();

		for (int I = 0; I <= m.getAtomCount() - 1; I++) { // start of atom 1-
			// need to
			// start process with each
			// atom to get all possible
			// paths

			List ca = m.getConnectedAtomsList(m.getAtom(I));

			if (ca.size() == M) {
				LinkedList al = new LinkedList();

				al.add(new Integer(m.getAtomNumber(m.getAtom(I))));
				for (int j = 0; j <= M - 1; j++) {					
					al.add(new Integer(m.getAtomNumber(((IAtom)ca.get(j)))));
				}
				MasterList.add(al);
			} else if (ca.size() > M) { 

				LinkedList mlist=FindAllPossibleCombos(M,ca.size());
							
				for (int i=0;i<=mlist.size()-1;i++) {
					LinkedList al = new LinkedList();
					al.add(new Integer(m.getAtomNumber(m.getAtom(I))));
					List l=Utilities.Parse((String)mlist.get(i),"\t");
					for (int j=0;j<=l.size()-1;j++) {
						String snum=(String)l.get(j);
						int num=Integer.parseInt(snum);						
						al.add(new Integer(m.getAtomNumber(((IAtom)ca.get(num)))));
					}
					MasterList.add(al);
					
				}
				
								
				/*if (M==3 && BoundAtoms.length==4) {
					// do all possibilities:
					LinkedList al1 = new LinkedList();
					al1.add(m.getAtomNumber(m.getAtom(I)) + "");
					al1.add(m.getAtomNumber(BoundAtoms[0]) + "");
					al1.add(m.getAtomNumber(BoundAtoms[1]) + "");
					al1.add(m.getAtomNumber(BoundAtoms[2]) + "");
					
					LinkedList al2 = new LinkedList();
					al2.add(m.getAtomNumber(m.getAtom(I)) + "");
					al2.add(m.getAtomNumber(BoundAtoms[0]) + "");
					al2.add(m.getAtomNumber(BoundAtoms[2]) + "");
					al2.add(m.getAtomNumber(BoundAtoms[3]) + "");
					
					LinkedList al3 = new LinkedList();
					al3.add(m.getAtomNumber(m.getAtom(I)) + "");
					al3.add(m.getAtomNumber(BoundAtoms[1]) + "");
					al3.add(m.getAtomNumber(BoundAtoms[2]) + "");
					al3.add(m.getAtomNumber(BoundAtoms[3]) + "");
					
					LinkedList al4 = new LinkedList();
					al4.add(m.getAtomNumber(m.getAtom(I)) + "");
					al4.add(m.getAtomNumber(BoundAtoms[0]) + "");
					al4.add(m.getAtomNumber(BoundAtoms[1]) + "");
					al4.add(m.getAtomNumber(BoundAtoms[3]) + "");
					
					MasterList.add(al1);
					MasterList.add(al2);
					MasterList.add(al3);
					MasterList.add(al4);
					
					
				} else if (M==3 && BoundAtoms.length==5) {
					//TODO
				} else if (M==4 && BoundAtoms.length==5) {
					//TODO
				} else if (M==3 && BoundAtoms.length==6) {
//					TODO
				} else if (M==4 && BoundAtoms.length==6) {
//					TODO
				}
*/
			}

		} // end i atom for loop

		// PrintList(MasterList);

		
		/*
		
		// duplicate checking is not necessary for clusters so following was commented out
		 
		FlagDuplicates(MasterList, 2);
		
		LinkedList MasterList2=new LinkedList();
		
		for (int i=0;i<=MasterList.size()-1;i++) {
			LinkedList AL=(LinkedList)MasterList.get(i);
			
			if (!AL.get(0).equals("duplicate")) {
				MasterList2.add(AL);
			}
			
		}
		
		if (MasterList.size()!=MasterList2.size()) System.out.println("found cluster duplicate");
		//System.out.println("Clusters:"+MasterList.size()+"\t"+MasterList2.size());
		
						
		return MasterList2;
		*/
		
		
		return MasterList;

	}
	
	static LinkedList FindAllPossibleCombos(int M,int NBound) {
		String s1="",s2="",s3="",s4="";
		
		LinkedList mlist=new LinkedList();
						
		for (int i=0;i<=NBound-1;i++) {
			
			s1=i+"";
			
			for (int j=i+1;j<=NBound-1;j++) {
								
				s2=s1+"\t"+j;
								
				for (int k=j+1;k<=NBound-1;k++) {
										
					s3=s2+"\t"+k;
					
					if (M==3) {
						mlist.add(s3);
					} else if (M==4) {
						
						for (int l=k+1;l<=NBound-1;l++) {
							s4=s3+"\t"+l;
							mlist.add(s4);
						} // end l for loop
						
					}// end M==4 else if
					
					
				} // end k for loop
				
			} // end j for loop
			
		} // end i for loop
		
		return mlist;
		
	}
	
		
	
	public static void FlagDuplicates(LinkedList MasterList, int bob) {

		
		Iloop: for (int i = 0; i <= MasterList.size() - 1; i++) {
			LinkedList al1 = (LinkedList) MasterList.get(i);

			// changed lower bound of j from 0 to i+1 to speed it up
			for (int j = i + 1; j <= MasterList.size() - 1; j++) {

				LinkedList al2 = (LinkedList) MasterList.get(j);
				boolean duplicate = false;

				if (bob == 1)
					duplicate = CompareLists(al1, al2);
				else if (bob == 2)
					duplicate = CompareLists2(al1, al2);

				if (duplicate) {
					al1.set(0,"duplicate");
					
					//MasterList.remove(i);
					//i--;
					continue Iloop;
				}

			}

		}

	}
	
	public static void RemoveDuplicates(LinkedList MasterList, int bob) {
		
		// essentially each path has 1 duplicate since the paths are the reverse
		// of each other when we start at the atom at the end of path

		// create a new MasterList with strings instead of Lists for each element
		LinkedList MasterList2=new LinkedList();
		
		
		ListIterator listIterator=MasterList.listIterator();
		
		// convert master list (list of lists) to list of strings:
		
		while (listIterator.hasNext()) {
			
			LinkedList ll=(LinkedList)listIterator.next();
			
			String strFirst=(String)ll.get(0);
			int First=Integer.parseInt(strFirst);
			
			String strLast=(String)ll.get(ll.size()-1);
			int Last=Integer.parseInt(strLast);
			
			String strPath="";
						
			
			if (First<Last) {
				for (int i=0;i<=ll.size()-2;i++) {
					strPath+=ll.get(i)+"-";
				}
				strPath+=ll.get(ll.size()-1);
			} else {
				for (int i=ll.size()-1;i>=1;i--) {
					strPath+=ll.get(i)+"-";
				}
				strPath+=ll.get(0);								
			}
			
			MasterList2.add(strPath);
			
		}
		
		// sort the list of strings:
		Collections.sort(MasterList2);
		
		ListIterator li2=MasterList2.listIterator();
		
		//System.out.println(MasterList.size());
		
	
		// remove the duplicates:
		while (li2.hasNext()) {
			String p1=(String)li2.next();
			
			//System.out.println("p1="+p1);
						
			if (li2.hasNext()) {
				
				String p2=(String)li2.next();
				//System.out.println("p2="+p2);
				
				if (p1.equals(p2)) {
					//MasterList2.remove(); // messes up list iterator!
					li2.remove();											
				}			
			}
			
		}
		// clear out masterlist
		MasterList.clear();
		
		li2=MasterList2.listIterator();
		
		// convert master list back to list of linkedlists:
		
		while (li2.hasNext()) {
			LinkedList ll=Utilities.Parse((String)li2.next(),"-");
			MasterList.add(ll);			
		}
				
		
		
		/*
		if (true) return;
		
		// System.out.println(MasterList.size());

		
		//ListIterator LI=MasterList.listIterator();
				
		Iloop: 									
			//while (LI.hasNext()) {
			
			for (int i = 0; i <= MasterList.size() - 1; i++) {
				LinkedList al1 = (LinkedList) MasterList.get(i);
				
				ListIterator li=MasterList.listIterator(i+1);
				
				while (li.hasNext()) {
				// changed lower bound of j from 0 to i+1 to speed it up
				//for (int j = i + 1; j <= MasterList.size() - 1; j++) {
					
					//if (i==0) System.out.println(li.nextIndex());
					
					LinkedList al2 = (LinkedList) li.next();
					boolean duplicate = false;
					
					if (bob == 1)
						duplicate = CompareLists(al1, al2);
					else if (bob == 2)
						duplicate = CompareLists2(al1, al2);
					
					if (duplicate) {						
						li.remove();							
						continue Iloop;
					}
					
				}
				
			}
*/
	}
	
	
	public static boolean CompareLists(LinkedList<Integer> al1, LinkedList<Integer> al2) {

		// two paths are equal if order of atoms is the same in either
		// forward or reverse directions

		boolean duplicate = true;

		LinkedList BondList1 = new LinkedList();
		LinkedList BondList2 = new LinkedList();

		int first1 = al1.get(0);
		int first2 = al2.get(0);
		int last2 = al2.get(al2.size() - 1);

		if (first1 != first2 && first1 != last2) {// if first element of path
			// 1 doesnt equal first or
			// last element of path 2
			// stop we dont have
			// duplicate
			duplicate = false;
			return duplicate;
		}

		String bob="";
		
		
		for (int i = 0; i <= al1.size() - 1; i++) {
			int one = al1.get(i);
			int two;

			if (first1 == first2) {
				two = al2.get(i);
			} else { // path2 might be reverse of path 1
				two = al2.get(al2.size() - 1 - i);
			}

			if (one != two) {
				duplicate = false;
				return duplicate;
			}

		}

		/*
		 * for (int i=0;i<=al1.size()-2;i++) { int
		 * one=Integer.parseInt((String)al1.get(i)); int
		 * two=Integer.parseInt((String)al1.get(i+1));
		 * 
		 * int three=Integer.parseInt((String)al2.get(i)); int
		 * four=Integer.parseInt((String)al2.get(i+1));
		 * 
		 * if (one < two) { BondList1.add(one+"-"+two); } else {
		 * BondList1.add(two+"-"+one); }
		 * 
		 * if (three < four) { BondList2.add(three+"-"+four); } else {
		 * BondList2.add(four+"-"+three); }
		 *  } // end i for list
		 * 
		 * java.util.Collections.sort(BondList1);
		 * java.util.Collections.sort(BondList2);
		 * 
		 * 
		 * 
		 * for (int i=0;i<=BondList1.size()-1;i++) { String
		 * one=(String)BondList1.get(i); String two=(String)BondList2.get(i);
		 * 
		 * if (!one.equals(two)) { duplicate=false; return duplicate; } }
		 */

		return duplicate;
	}

	public static boolean CompareLists2(LinkedList al1, LinkedList al2) {

		// two path-clusters are equal if all bonds are the same
		// al1 and al2 are LinkedLists of bond strings

		boolean duplicate = true;

		LinkedList al3 = (LinkedList) al1.clone();
		LinkedList al4 = (LinkedList) al2.clone();

		java.util.Collections.sort(al3);
		java.util.Collections.sort(al4);

		// LinkedList BondList2=new LinkedList();

		for (int i = 0; i <= al3.size() - 1; i++) {
			String one = (String) al3.get(i);
			String two = (String) al4.get(i);
			if (!one.equals(two)) {
				duplicate = false;
				// System.out.println(duplicate);
				return duplicate;
			}
		}

		return duplicate;
	}

	public static LinkedList FindPathClusters(IAtomContainer m) {

		LinkedList MasterList = new LinkedList();
		LinkedList MasterBondList = new LinkedList();

		// first find all third order clusters

		for (int I = 0; I <= m.getAtomCount() - 1; I++) { // start of atom 1-
			// need to
			// start process with each
			// atom to get all possible
			// paths

			List ca = m.getConnectedAtomsList(m.getAtom(I));

			String strBond;

			if (ca.size() == 3) {

				LinkedList al = new LinkedList();
				LinkedList BondList = new LinkedList();

				int ii = m.getAtomNumber(m.getAtom(I));

				al.add(new Integer(ii));
				for (int j = 0; j <= 3 - 1; j++) {
					int jj = m.getAtomNumber((IAtom)ca.get(j));
					al.add(new Integer(jj));

					if (ii < jj) {
						strBond = ii + "-" + jj;
					} else {
						strBond = jj + "-" + ii;
					}
					BondList.add(strBond);

				}

				MasterBondList.add(BondList);
				MasterList.add(al);

			} else if (ca.size() >= 4) { // case where finding third
				// order clusters and have more than 3 atoms attached
			
				int ii = m.getAtomNumber(m.getAtom(I));

				List mlist=FindAllPossibleCombos(3,ca.size());
				
				for (int i=0;i<=mlist.size()-1;i++) {
					LinkedList al = new LinkedList();
					LinkedList BondList = new LinkedList();
					
					al.add(new Integer(ii));
					List l=Utilities.Parse((String)mlist.get(i),"\t");
					for (int j=0;j<=l.size()-1;j++) {
						String snum=(String)l.get(j);
						int num=Integer.parseInt(snum);						
						int jj=m.getAtomNumber((IAtom)ca.get(num));
						al.add(new Integer(jj));
						
						if (ii < jj) {
							strBond = ii + "-" + jj;
						} else {
							strBond = jj + "-" + ii;
						}
						BondList.add(strBond);
						
					}
					MasterList.add(al);
					MasterBondList.add(BondList);
					
				}
											
				/*				
				int[][] BobArray = { { 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 },
						{ 1, 2, 3 } };

				for (int i = 0; i <= 3; i++) {
					LinkedList al1 = new LinkedList();
					LinkedList BondList = new LinkedList();

					al1.add(ii + "");

					for (int j = 0; j <= 2; j++) {

						// System.out.println(BobArray[i][j]);
						int jj = m.getAtomNumber(BoundAtoms[BobArray[i][j]]);
						al1.add(jj + "");
						if (ii < jj) {
							strBond = ii + "-" + jj;
						} else {
							strBond = jj + "-" + ii;
						}
						BondList.add(strBond);

					}

					MasterList.add(al1);
					MasterBondList.add(BondList);
				}
				
				*/

			} // end else if (BoundAtoms.length >= 4)

		} // end I for loop

		// PrintList(MasterList);
		// PrintList(MasterBondList);

		// System.out.println("\n-----\n");

		int Size = MasterList.size() - 1;

		// now find all path clusters which are bound to each of the third order
		// path clusters:
		for (int i = 0; i <= Size; i++) {
			LinkedList <Integer>al = (LinkedList) MasterList.get(i);
			LinkedList BondList = (LinkedList) MasterBondList.get(i);

			int atomnum0 = al.get(0);
			for (int j = 1; j <= al.size() - 1; j++) {
				int atomnum = al.get(j);

				List BoundAtoms = m.getConnectedAtomsList(m.getAtom(atomnum));

				for (int k = 0; k <= BoundAtoms.size() - 1; k++) {

					if (m.getAtomNumber((IAtom)BoundAtoms.get(k)) != atomnum0) {

						LinkedList al2 = (LinkedList) al.clone();
						LinkedList BondList2 = (LinkedList) BondList.clone();
						int ii = atomnum;
						int jj = m.getAtomNumber((IAtom)BoundAtoms.get(k));

						al2.add(new Integer(jj));
						MasterList.add(al2);

						String strBond;
						if (ii < jj) {
							strBond = ii + "-" + jj;
						} else {
							strBond = jj + "-" + ii;
						}
						BondList2.add(strBond);
						MasterBondList.add(BondList2);

					}
				}

			}

		}

		for (int i = 0; i <= Size; i++) {
			MasterList.remove(0);
			MasterBondList.remove(0);
		}

		// PrintList(MasterBondList);

		// check to see if any lists have same set of bonds:

		Iloop: for (int i = 0; i <= MasterBondList.size() - 1; i++) {
			LinkedList al1 = (LinkedList) MasterBondList.get(i);

			for (int j = 0; j <= MasterBondList.size() - 1; j++) {
				if (i != j) {
					LinkedList al2 = (LinkedList) MasterBondList.get(j);
					boolean duplicate;

					duplicate = CompareLists2(al1, al2);

					if (duplicate) {
						MasterList.remove(i);
						MasterBondList.remove(i);
						i--;
						continue Iloop;
					}
				}
			}

		}

		// check to see if any atoms are repeated in each list:

		for (int i = 0; i <= MasterList.size() - 1; i++) {

			LinkedList myList = (LinkedList) MasterList.get(i);
			boolean repeat = CheckForRepeatedAtoms(myList);

			if (repeat) {
				MasterList.remove(i);
				i--;
			}
		}

		// PrintList(MasterList);

		return MasterList;

	}
	
	public static boolean CheckForRepeatedAtoms(LinkedList myList) {

		boolean repeat = false;

		// make a copy of myList:
		LinkedList <Integer>bob = (LinkedList) myList.clone();
		java.util.Collections.sort(bob);

		for (int i = 0; i <= bob.size() - 2; i++) {
			int one = bob.get(i);
			int two = bob.get(i + 1);
			if (one == two) {
				repeat = true;

				return repeat;
			}

		}

		return repeat;

	}

	/**
	 * returns MasterList which is an array of linkedlists each element of the
	 * array represents the paths for a given path length each path is given by
	 * the atom numbers of the atoms in the path (not an array of atoms) Now
	 * atom numbers are stored as Integer instead of String to save memory
	 * 
	 * @param M
	 * @param m
	 * @return
	 */
	public static LinkedList[] FindPaths4(int M,IAtomContainer m) {
		Runtime s_runtime = Runtime.getRuntime ();
		long start = System.currentTimeMillis();
		int ArraySize;

		if (M < m.getAtomCount())
			ArraySize = m.getAtomCount() + 1;
		else
			ArraySize = M + 2;

		LinkedList[] MasterList = new LinkedList[ArraySize];

		for (int i = 1; i <= ArraySize - 1; i++) {
			MasterList[i] = new LinkedList();
		}
		
//		ToxPredictor.MyDescriptors.PathTools pt=new ToxPredictor.MyDescriptors.PathTools();
		
//		 unique paths
        for (int i = 0; i < m.getAtomCount()- 1; i++) {
            IAtom a = m.getAtom(i);
            for (int j = i + 1; j < m.getAtomCount(); j++) {
                IAtom b = m.getAtom(j);
                
//                List list=PathTools.getAllPaths(m, a, b);
//                List list=pt.getAllPaths(m, a, b);
                List list=org.openscience.cdk.graph.PathTools.getAllPaths(m, a, b);
                
                
                if ((System.currentTimeMillis()-start) > TESTConstants.pathGenerationTimeout) {
                    throw new RuntimeException(
                            String.format("Timeout %s ms while generating paths for %s.", 
                                    TESTConstants.pathGenerationTimeout, 
                                    (String)m.getProperty("CAS")));
                }
                
                for (int k=0;k<list.size();k++) {
                	List list2=(List)list.get(k);
//                	MasterList[list2.size()-1].add(list2);

                	LinkedList list3=new LinkedList();
                	for (int l=0;l<list2.size();l++) {
                		IAtom ia=(IAtom)list2.get(l);
                		list3.add(new Integer(m.getAtomNumber(ia)));
                		
//                		System.out.println(s_runtime.freeMemory());
                		if (s_runtime.freeMemory()<10000) {
//                			System.out.println(s_runtime.freeMemory());
                			return null;
                		}
                		
                	}
                	
                	MasterList[list2.size()-1].add(list3);
                	
                }
                
                
            }
        }
        
        

//        System.out.println("\n\n"+m.getProperty("CAS"));
//        
//        for (int i=1;i<MasterList.length;i++) {
//			
//			List list=MasterList[i];
//
//			System.out.println("i="+i);
//			for (int j=0;j<list.size();j++) {
//				List list2=(List)list.get(j);
//				
//				for (int k=0;k<list2.size();k++) {					
//					System.out.print(list2.get(k)+"\t");
//				}
//				System.out.print("\n");
//			}
//		}
		
		
		return MasterList;
	}
	
	
	public static LinkedList[] FindPaths3(int M,IAtomContainer m) {

		// returns MasterList which is an array of linkedlists
		// each element of the array represents the paths for a given path length
		// each path is given by the atom numbers of the atoms in the path (not an array of atoms)		
		//old method
		// this method stores atom numbers as strings- which is incompatible with rest of descriptor calculations

		try {

			// LinkedList MasterList=new LinkedList();

			//	 list for given starting  atom I:
			LinkedList PathList = new LinkedList(); 
			
			
			int ArraySize;

			if (M < m.getAtomCount())
				ArraySize = m.getAtomCount() + 1;
			else
				ArraySize = M + 2;

			LinkedList[] MasterList = new LinkedList[ArraySize];

			for (int i = 1; i <= ArraySize - 1; i++) {
				MasterList[i] = new LinkedList();
			}

			// for (int I=1;I<=m.NATOMS;I++) { // start of atom 1- need to start
			// process with each atom to get all possible paths
			for (int I = 0; I <= m.getAtomCount() - 1; I++) { // start of atom
				// 1- need
				// to start process with
				// each atom to get all
				// possible paths

				// System.out.println("start atom "+(I+1)+" of
				// "+m.getAtomCount());

				PathList.clear();

				LinkedList newList = new LinkedList();
				newList.add(I + "");
				PathList.add(newList);

				LinkedList l = null;

				while (true) {
					// for (int Z=1;Z<=2;Z++) {

					int CurrentPosition = PathList.size() - 1;

					KLoop: for (int K = 0; K <= CurrentPosition; K++) {

						LinkedList myList = (LinkedList) PathList.get(K);

						int CurrentAtomNum = Integer.parseInt((String) myList
								.get(myList.size() - 1));
						int PreviousAtomNum;

						if (myList.size() - 2 >= 0) {
							PreviousAtomNum = Integer.parseInt((String) myList
									.get(myList.size() - 2));
						} else {
							PreviousAtomNum = -1;
						}

						// now find atoms bound current atom:
						List BoundAtoms = m.getConnectedAtomsList(m
								.getAtom(CurrentAtomNum));

						// need to make sure dont add previous atom

						if (BoundAtoms.size() > 0) {
							for (int i = 0; i < BoundAtoms.size(); i++) {
								LinkedList templist = (LinkedList) (PathList
										.get(K));
								LinkedList newList2 = (LinkedList) templist
										.clone();
								IAtom a = (IAtom)BoundAtoms.get(i);

								if (m.getAtomNumber(a) != PreviousAtomNum) {
									newList2.add(m.getAtomNumber(a) + "");
									PathList.add(newList2);

								}

							} // end J loop

						}

					} // end K for loop

					// PrintList(PathList);

					// delete previous list: (paths with one less atom in list)
					for (int i = 0; i <= CurrentPosition; i++) {
						PathList.remove(0);
					}

					// PrintList(PathList);

					if (PathList.size() == 0)
						break;

					l = (LinkedList) PathList.get(0);

					// PrintList(PathList);

					// check to see if any atoms are repeated in each list:
					for (int i = 0; i <= PathList.size() - 1; i++) {

						LinkedList myList = (LinkedList) PathList.get(i);
						boolean repeat = CheckForRepeatedAtoms(myList);

						if (repeat) {
							PathList.remove(i);
							i--;
						}
					}

					// System.out.println(PathList.size());
					// PrintList(PathList);

					for (int i = 0; i <= PathList.size() - 1; i++) {
						LinkedList al = (LinkedList) PathList.get(i);
						LinkedList alclone = (LinkedList) al.clone();
						int size = alclone.size();
						MasterList[size - 1].add(alclone);
						// System.out.println(size-1+"\t"+MasterList[size-1].size());
					}

					// PrintList(PathList);

					if (PathList.size() > 0) {
						l = (LinkedList) PathList.get(0);
						// System.out.println(l.size());

						if (l.size() == m.getAtomCount()) {
							break;
						}

					} else {
						break;
					}

				} // end while true;

			} // end I loop for loop

			// System.out.println("end I for loop");
			// PrintList(MasterList);

			// now need to check master list for duplicates:

			// checking for path duplicates:
			// System.out.print("Checking for path duplicates...");

			
			double time1=System.currentTimeMillis()/1000.0;
			
			
			
			
			//System.out.println("");
			
			for (int i = 1; i <= ArraySize-1; i++) {			
			//for (int i = 1; i <= M; i++) { // for speed stop at M
				
				 //System.out.println(i+" before:"+MasterList[i].size());
				//FlagDuplicates(MasterList[i], 1); // note: flagging then only adding non duplicates had only a very small increase in performance (1% for strychnine)
				RemoveDuplicates(MasterList[i],1);
				 //System.out.println(i+" after:"+MasterList[i].size());
			}
			
			
			/*
			// create a new master list and only add non duplicates
			LinkedList [] MasterList2=new LinkedList [MasterList.length];
									
			for (int i=1;i<=MasterList2.length-1;i++) {
				LinkedList al=new LinkedList();
				MasterList2[i]=al;
				
				LinkedList AL=(LinkedList)MasterList[i];
				//System.out.println(i);
				for (int j=0;j<=AL.size()-1;j++) {
					LinkedList AL2=(LinkedList)AL.get(j);
					if (!AL2.get(0).equals("duplicate")) {
						MasterList2[i].add(AL2);
					}
					
				}
				
			}
			*/
			
			double time2=System.currentTimeMillis()/1000.0;
			
			//System.out.println("duplicate checking took "+(time2-time1)+" seconds");
			
			//for (int i=1;i<=MasterList.length-1;i++) {
//				System.out.println(i+"\t"+MasterList[i].size()+"\t"+MasterList2[i].size());
			//}
			
			// System.out.print("done\n");

			// PrintList(MasterList);

			//return MasterList2;
			return MasterList;

		} catch (Exception e) {
			e.printStackTrace();
			return null;

		}

	}
	
	
	/*
	 * ArrayList FindChains(int M,int [][]IC) { try {
	 * 
	 * ArrayList MasterList=new ArrayList(); ArrayList PathList=new ArrayList(); //
	 * list for given starting atom I
	 * 
	 * if (m.NATOMS<M+1) return MasterList;
	 * 
	 * 
	 * //for (int I=1;I<=m.NATOMS;I++) { // start of atom 1- need to start
	 * process with each atom to get all possible paths for (int I=1;I<=m.NATOMS;I++) { //
	 * start of atom 1- need to start process with each atom to get all possible
	 * paths PathList.clear();
	 * 
	 * 
	 * ArrayList newList=new ArrayList(); newList.add(I+"");
	 * PathList.add(newList);
	 * 
	 * ArrayList l=null;
	 * 
	 * 
	 * while (true) { //for (int Z=1;Z<=2;Z++) {
	 * 
	 * int CurrentPosition=PathList.size()-1;
	 * 
	 * KLoop: for (int K=0;K<=CurrentPosition;K++) {
	 * 
	 * ArrayList myList=(ArrayList)PathList.get(K);
	 * 
	 * int CurrentAtomNum=Integer.parseInt((String)myList.get(myList.size()-1));
	 * int PreviousAtomNum;
	 * 
	 * if (myList.size()-2>=0) {
	 * PreviousAtomNum=Integer.parseInt((String)myList.get(myList.size()-2)); }
	 * else { PreviousAtomNum=-1; } //System.out.println(K+"\t"+CurrentAtomNum); //
	 * now find atoms bound current atom: ArrayList BoundAtoms=new ArrayList();
	 * 
	 * for (int J=1;J<=m.NATOMS;J++) { if (J!=CurrentAtomNum) { if
	 * (IC[CurrentAtomNum][J]!=0) { if (J!=PreviousAtomNum) { // dont go back to
	 * previous atom BoundAtoms.add(J+""); } } } } // end J loop
	 * 
	 * if (BoundAtoms.size()>0) { for (int i=0;i<=BoundAtoms.size()-1;i++) {
	 * ArrayList templist=(ArrayList)(PathList.get(K)); ArrayList
	 * newList2=(ArrayList)templist.clone(); newList2.add(BoundAtoms.get(i));
	 * PathList.add(newList2); } // end J loop } } // end K for loop
	 * 
	 * 
	 * //PrintList(PathList); // delete previous list: (paths with one less atom
	 * in list) for (int i=0;i<=CurrentPosition;i++) { PathList.remove(0); }
	 * 
	 * //PrintList(PathList);
	 * 
	 * if (PathList.size()==0) break;
	 * 
	 * l=(ArrayList)PathList.get(0);
	 * 
	 * //PrintList(PathList); // now need to check to see if have any repeats
	 * other than first and last atom: iloop: for (int i=0;i<=PathList.size()-1;i++) {
	 * 
	 * ArrayList myList=(ArrayList)PathList.get(i);
	 * 
	 * for (int j=0;j<=myList.size()-1;j++) {
	 * 
	 * int one=Integer.parseInt((String)myList.get(j));
	 * 
	 * for (int k=0;k<=myList.size()-1;k++) { int
	 * two=Integer.parseInt((String)myList.get(k));
	 * 
	 * if (one==two) { if (j!=k && k<myList.size()-1) { PathList.remove(i);
	 * i=-1; continue iloop; }
	 * 
	 * if (j==0 && k==myList.size()-1 && myList.size()==M+1) { ArrayList
	 * myList2=(ArrayList)myList.clone(); MasterList.add(myList2);
	 * PathList.remove(i); i=-1; continue iloop; } } // end if one=two } // end
	 * k for loop } // end j for loop } // end i for loop
	 * 
	 * 
	 * //PrintList(PathList);
	 * 
	 * if (PathList.size()>0) { l=(ArrayList)PathList.get(0);
	 * //System.out.println(l.size());
	 * 
	 * if (l.size()==M+1) { break; } } else { break; } } // end while true; } //
	 * end I loop for loop
	 * 
	 * 
	 * //PrintList(MasterList); // now need to check master list for duplicates:
	 * this.RemoveDuplicates(MasterList); // now trim off last atom off each
	 * list (since its a duplicate): for (int j=0;j<=MasterList.size()-1;j++) {
	 * ArrayList al=(ArrayList)MasterList.get(j); al.remove(al.size()-1); }
	 * 
	 * //PrintList(MasterList);
	 * 
	 * 
	 * return MasterList; } catch (Exception ex) { ex.printStackTrace(); return
	 * new ArrayList(); } }
	 * 
	 * 
	 */

	
	
	public static LinkedList[] FindFragmentPaths(IAtomContainer m,int atomnum, String[] Fragment,
			LinkedList[] paths) {

		int ArraySize;
		int M = 10;

		if (M < m.getAtomCount())
			ArraySize = m.getAtomCount() + 1;
		else
			ArraySize = M + 2;

		LinkedList[] MasterList = new LinkedList[ArraySize];

		// first need to find all fragment paths that start with the current
		// atom

		for (int i = 1; i <= ArraySize - 1; i++) {
			MasterList[i] = new LinkedList();

			LinkedList al = (LinkedList) paths[i];

			for (int j = 0; j <= al.size() - 1; j++) {

				LinkedList al2 = (LinkedList) al.get(j);

				int first = Integer.parseInt((String) al2.get(0));
				int last = Integer.parseInt((String) al2.get(al2.size() - 1));

				if (atomnum == first || atomnum == last) {

					String fragmentpath = "";
					for (int k = 1; k <= al2.size() - 1; k++) {
						// start at one since dont want to include start atom
						
						int atomnum2;
						if (atomnum == first)
							atomnum2 = Integer.parseInt((String) al2.get(k));
						else
							atomnum2 = Integer.parseInt((String) al2.get(al2
									.size()
									- 1 - k));

						
						String f=Fragment[atomnum2];
												
						if (f instanceof String) {
							// replace since topologically equivalent:
							if (f.equals("SsOm"))
								f = "SdO";
						}

						if (k < al2.size() - 1)
							fragmentpath += f + "-";
						else
							fragmentpath += f;
					}
					MasterList[i].add(fragmentpath);
				}

			} // end j for loop

			Lloop: for (int l = 0; l <= MasterList[i].size() - 1; l++) {
				String one = (String) MasterList[i].get(l);

				// changed lower bound of j from 0 to i+1 to speed it up
				for (int mm = l + 1; mm <= MasterList[i].size() - 1; mm++) {
					String two = (String) MasterList[i].get(mm);

					boolean duplicate = false;

					if (one.equals(two)) {
						MasterList[i].remove(l);
						l--;
						continue Lloop;
					}

				}

			}

			// sort lists:

			LinkedList bob = (LinkedList) MasterList[i];
			Collections.sort(bob, Collator.getInstance()); // need collator to
			// sort strings with
			// letters properly
			/*
			 * for (int z=0;z<=MasterList[i].size()-1;z++) {
			 * System.out.println(MasterList[i].get(z)); }
			 */

		} // end i for loop

		return MasterList;
	}

	public static boolean CompareFragmentPaths(IAtomContainer m,LinkedList[] fpaths1, LinkedList[] fpaths2) {

		// System.out.println("enter compare frag");
		boolean same = true;

		int ArraySize;
		int M = 10;

		if (M < m.getAtomCount())
			ArraySize = m.getAtomCount() + 1;
		else
			ArraySize = M + 2;

		// first need to find all fragment paths that start with the current
		// atom

		for (int i = 1; i <= ArraySize - 1; i++) {
			LinkedList al1 = fpaths1[i];
			LinkedList al2 = fpaths2[i];

			if (al1.size() != al2.size()) {
				same = false;
				return same;
			}

			for (int j = 0; j <= al1.size() - 1; j++) {
				String s1 = (String) al1.get(j);
				String s2 = (String) al2.get(j);

				if (!s1.equals(s2)) {

					same = false;
					return same;
				}
			}

		}

		return same;

	}

	
	public static void PrintList(LinkedList l) {

		if (l.size() == 0)
			return;

		LinkedList l2 = (LinkedList) l.get(0);
		// System.out.println(l2.size()+":");

		for (int i = 0; i <= l.size() - 1; i++) {
			LinkedList myList = (LinkedList) l.get(i);

			for (int j = 0; j <= myList.size() - 1; j++) {
				System.out.print(myList.get(j) + "\t");
			}
			System.out.print("\n");
		}
		System.out.print("\n");

	}

	public static void PrintList2(LinkedList l) {

		if (l.size() == 0)
			return;

		for (int j = 0; j <= l.size() - 1; j++) {
			System.out.print(l.get(j) + "\t");
		}
		System.out.print("\n");

	}

	public static void main(String[] args) {
		IAtomContainer m=null;
	     SmilesParser   sp  = new SmilesParser(DefaultChemObjectBuilder.getInstance());

		
		String smiles="c1ccccc1";
		
		
		smiles="O=C(N(c(c(C1(C(N(C2)CC(C3C4C5OC6)=C6)C3)C2)ccc7)c7)C14)C5";
		try {
		     m   = (IAtomContainer)sp.parseSmiles(smiles);

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println("here");
		double time1=System.currentTimeMillis()/1000.0;		
		PathFinder.FindPaths3(10,m);
		double time2=System.currentTimeMillis()/1000.0;
		System.out.println(time2-time1);
	}
	
	/*
	 * private void CalculateDistanceMatrix(int[][] Distance, ArrayList[]
	 * PathList) { // this sub is now unnecessary- can calc using: // calculate
	 * adjacency matrix: //int[][] admat =
	 * org.openscience.cdk.graph.matrix.AdjacencyMatrix.getMatrix(m); //Distance =
	 * org.openscience.cdk.graph.PathTools.computeFloydAPSP(admat);
	 * 
	 * 
	 * for (int i = 0; i <= m.getAtomCount() - 1; i++) { for (int j = 0; j <=
	 * m.getAtomCount() - 1; j++) { Distance[i][j] = -1; } }
	 * 
	 * iloop: for (int i = 0; i <= m.getAtomCount() - 1; i++) { //
	 * System.out.println("i="+i);
	 * 
	 * jloop: for (int j = 0; j <= m.getAtomCount() - 1; j++) { //
	 * System.out.println("j="+j);
	 * 
	 * if (i == j) { Distance[i][j] = 0; Distance[j][i] = 0; continue jloop; }
	 * 
	 * if (Distance[i][j] > -1) continue jloop; //
	 * System.out.println(i+"\t"+j+"\t"+"here");
	 * 
	 * for (int pathlength = 1; pathlength <= m.getAtomCount(); pathlength++) {
	 * 
	 * ArrayList al = (ArrayList) PathList[pathlength];
	 * 
	 * for (int k = 0; k <= al.size() - 1; k++) {
	 * 
	 * ArrayList al2 = (ArrayList) al.get(k);
	 * 
	 * if (al2.size() > 0) { int ii = Integer.parseInt((String) (al2.get(0)));
	 * int jj = Integer.parseInt((String) (al2.get(al2 .size() - 1)));
	 * 
	 * if ((i == ii && j == jj) || (i == jj && j == ii)) { Distance[i][j] =
	 * al2.size() - 1; Distance[j][i] = Distance[i][j]; //
	 * System.out.println(i+"\t"+j+"\t"+Distance[i][j]); continue jloop; } } }
	 * 
	 * }// end path length loop } // end j for loop }// end i for loop // for
	 * (int i=0;i<=m.getAtomCount()-1;i++) { for (int j=0;j<=m.getAtomCount()-1;j++) {
	 * //System.out.print(Distance[i][j]+"\t"); } System.out.print("\n"); } }
	 */
} // end class

