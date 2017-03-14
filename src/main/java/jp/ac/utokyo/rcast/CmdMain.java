package jp.ac.utokyo.rcast;

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


import java.io.File;
import java.io.IOException;
import java.util.Optional;
import jp.ac.utokyo.rcast.realign.Realignment;
import org.apache.commons.lang3.ArrayUtils;

public class CmdMain {
  public static void main(final String[] args) throws IOException, InterruptedException {
    if (args == null || args.length == 0 || args[0] == null) {
      printMessage();
      return;
    }

    switch (args[0]) {
      case "realign":
        Realignment.main(ArrayUtils.subarray(args, 1, args.length));
        break;
      default:
        printMessage();
        break;
    }
  }

  private static void printMessage() {
    final Optional<String> version = Optional.ofNullable(CmdMain.class.getPackage().getImplementationVersion());
    final Optional<String> title = Optional.ofNullable(CmdMain.class.getPackage().getImplementationTitle());
    final String jar = new File(CmdMain.class.getProtectionDomain().getCodeSource().getLocation().toString()).getName();
    System.err.println(title.flatMap(t->version.map(v->t+" version "+v)).orElse("dualRealign"));
    System.err.println("usage: java -jar " + jar + " {realign} <command> options");
  }
}
