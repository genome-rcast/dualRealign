/*
 * Copyright Hiroki Ueda
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


package srma.test;

import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.realign.Realignment;

public class TestRealgin {
	
	

	
	public static void main(String[] arg){
		

		
		String normalbamf = "/data/users/yamamoto/TodaiPanel/bam/PLC-TK-3N_TDv3_genome.bam";
		String tumorbamf = "/data/users/yamamoto/TodaiPanel/bam/PLC-TK-3TA_TDv3_genome.bam";
		
		
		String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		
		String targetRegion = "/data/users/yamamoto/TodaiPanel/target/S3035822_Covered.bed";
		
		
		String outdir = "/GLUSTER_DIST/data/users/ueda/toptest";
		
		List<String> l = new ArrayList<String>();
		
		add(l,"-n",normalbamf);
		add(l,"-t",tumorbamf);
		add(l,"-r",twobitref);
		
		add(l,"-ct",targetRegion);
	
		add(l,"-o",outdir);
		
		add(l,"-nt","8");
		
				
		
		String[] ar = l.toArray(new String[l.size()]);
		
		Realignment.main(ar);
		
	}

	private static void add(List<String> l, String s1, String s2) {
		l.add(s1);
		l.add(s2);
	}

}
