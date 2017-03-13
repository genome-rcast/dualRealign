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

import static org.junit.Assert.*;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.junit.Test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import jp.ac.utokyo.rcast.realign.Realignment;

public class CmdTestRealign {

  private static final String normal = "./src/test/resources/normal.bam";
  private static final String tumor = "./src/test/resources/tumor.bam";
  private static final String ref = "./src/test/resources/reference.2bit";
  private static final String target = "./src/test/resources/target.bed";

  @Test
  public void testRealign() throws Exception {
    // Arguments
    final Path tmpdir = Files.createTempDirectory("dualRealign");
    assertNotNull(tmpdir);

    // Execute
    final String[] argv = {"-n", normal, "-t", tumor, "-r", ref, "-ct", target, "-o", tmpdir.toString(), "-nt", "1"};
    Realignment.main(argv);

    // Checks
    final File n_result = Paths.get(tmpdir.toString(), "normal.bam_realign.bam").toFile();
    final File t_result = Paths.get(tmpdir.toString(), "tumor.bam_realign.bam").toFile();
    assertTrue(n_result.exists());
    assertTrue(t_result.exists());
    assertTrue(Paths.get(tmpdir.toString(), "normal.bam_realign.bai").toFile().exists());
    assertTrue(Paths.get(tmpdir.toString(), "tumor.bam_realign.bai").toFile().exists());

    final SamReaderFactory srf = SamReaderFactory.makeDefault();
    try(final SamReader n = srf.open(n_result);
        final SamReader t = srf.open(t_result);
        final SAMRecordIterator n_it = n.iterator();
        final SAMRecordIterator t_it = t.iterator()){
      assertEquals(3, toStream(n_it).count());
      assertEquals(3, toStream(t_it).count());
    }
  }

  private static <T> Stream<T> toStream(final Iterator<T> it) {
    final Iterable<T> iterable = () -> it;
    return StreamSupport.stream(iterable.spliterator(), false);
  }
}
