/*
 * Copyright Hiroki Ueda

 *  This program is free software; you can redistribute it and/or modify it under
 *	the terms of the GNU General Public License as published by the Free Software
 *	Foundation; either version 2 of the License, or (at your option) any later
 *	version.
	
 *	This program is distributed in the hope that it will be useful, but WITHOUT
 *	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *	FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 *	details.
	
 *	You should have received a copy of the GNU General Public License along with
 *	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 *	Place, Suite 330, Boston, MA 02111-1307 USA
 */


import java.io.IOException;


import jp.ac.utokyo.rcast.realign.Realignment;

public class CmdMain {

	//public static final String version = "Karkinos2.0 version 0.1 2013/4/22";
	public static final String version = "Karkinos2.0 version 0.3 2016/11/30";

	public static void main(String[] arg) throws Exception {

		//
		// arg = new String[]{"analysis"};
		if (arg == null || arg.length == 0) {
			printMsg();
		} else {

			String cmd = arg[0];
			String[] arg2 = getArg(arg);
			if (cmd.equals("realign")) {

				Realignment.main(arg2);	
				
			}else{

				printMsg();
			}

		}

	}

	private static String[] getArg(String[] arg) {
		String[] arg2 = new String[arg.length - 1];
		for (int n = 1; n < arg.length; n++) {
			arg2[n - 1] = arg[n];
			// System.out.println(arg[n]);
		}
		return arg2;
	}

	private static void printMsg() {
		System.out.println(version);
		System.out.println("usage: dualAlign.jar <command> options");
				System.out
		.println("		 realign	realignment normar,tumor bam file");
		
		

	}

}
