package umich.CNVfinder;

import htsjdk.samtools.*;
import sun.misc.GC;

import javax.management.Query;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.TreeMap;



public class Window {
    private Coord genomicInfo;
    private double GCcontent;
    private TreeMap<Long, Site> NTmap;
    private TreeMap<Long, Integer> RDmap;


    public Window(Coord _in) {
        genomicInfo = _in;
    }

    public TreeMap<Long, Integer> getRDmap() { return RDmap; }

    public double getPercentGC() {
        double ret = 0;
        double N = (double) NTmap.size();
        ret = GCcontent / N;
        return(ret);
    }

    public double getMedianReadCount() {
        double ret = 0;
        int[] rawReadDepth = new int[ RDmap.size() ];
        int i = 0;
        for(int rd : RDmap.values()) rawReadDepth[i++] = rd;

        if(rawReadDepth.length > 0) {

            Arrays.sort(rawReadDepth);
            int middle = rawReadDepth.length / 2;

            if(middle % 2 == 1) {
                ret = (double) rawReadDepth[middle];
            }
            else {
                double a = (double) (rawReadDepth[ (middle-1) ]);
                double b = (double) (rawReadDepth[ middle ]);
                ret = (a + b) / 2;
            }
        }
        return(ret);
    }


    // Function normalized the window read count by the GC content
    public void adjustReadCounts() {

        double m = Main.overall_median_RC;

        double bin = getBin(this.getPercentGC());

        double m_gc = Main.percentGCmap.get(bin);

        for(Long pos : RDmap.keySet()) {
            double ri = (double) RDmap.get(pos);

            double r_hat_i = ri * (m / m_gc );
            int new_r = (int) r_hat_i;
            RDmap.put(pos, new_r);
        }
    }


    private double getBin(double percentGC) {
        double ret = 0;
        double[] bins = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
        for(int j = 1; j < bins.length - 1; j++) {
            int i = j - 1;
            if( (percentGC >= bins[i]) && (percentGC < bins[j]) ) {
                ret = bins[i];
            }
        }
        return(ret);
    }


    // Prepare for reading the bam file.
    public void init() {

        NTmap = new TreeMap<>();

        for(long i = genomicInfo.start; i <= genomicInfo.end; i++) {
            Site curPos = new Site();
            NTmap.put(i, curPos);
        }
    }



    public void parseBAM(String bamFileName) throws IOException {

        File bamF = new File(bamFileName);

        if( !bamF.exists() ) {
           System.err.println("\nERROR: Unable to open " + bamFileName + "\n");
           System.exit(4);
        }

        // Open BAM file for reading
        final SamReader samReader = SamReaderFactory.makeDefault().open(bamF);
        int chrIndex = samReader.getFileHeader().getSequenceIndex(genomicInfo.chrom);
        if(chrIndex == -1) {
            String chrom = "chr" + genomicInfo.chrom;
            chrIndex = samReader.getFileHeader().getSequenceIndex(chrom);
        }

        final QueryInterval[] intervals = new QueryInterval[1];
        intervals[0] = new QueryInterval(chrIndex, (int) genomicInfo.start, (int) genomicInfo.end);

        // Iterate over the given BAM file
        // The 'b = true' argument requires the reads to be complete contained within the interval provided
        SAMRecordIterator iter = samReader.query(intervals, true);
        while(iter.hasNext()) {
            SAMRecord samRecord = iter.next();
            long start = (long) samRecord.getAlignmentStart();
            int MAPQ = samRecord.getMappingQuality();
            int flag = samRecord.getFlags();
            Cigar cigar = samRecord.getCigar();
            String rn = samRecord.getReadName().toUpperCase();

            // These are filters to apply to the read. If it fails any of them we don't record it
            if(MAPQ < Main.mapq_cutoff) continue;
            if( !(flag != 99) && !(flag != 147) ) continue;
            if( isSplitRead(cigar) ) continue;

            // Record the reads start and the GC content
            addReadDataByStart( rn, start, samRecord.getReadBases()[0] );
        }

        samReader.close();

        // Now that you are done parsing the BAM file, call all the sites as G/C
        recordGCcontent();
    }



    public void getReadDepthByWindow(int subWindowSize) {

        RDmap = new TreeMap<>();
        for(long i = genomicInfo.start; i <= (genomicInfo.end-subWindowSize); i += subWindowSize) {
            int RD = 0;
            for(long j = i; j < (i+subWindowSize); j++) {
                Site curSite = NTmap.get(j);
                RD += curSite.getReadDepth();
            }
            RDmap.put(i, RD);
        }

    }


    // Compute the GC content of this window.
    private void recordGCcontent() {
        GCcontent = 0;
        for(long i = genomicInfo.start; i <= genomicInfo.end; i++) {
            NTmap.get(i).callGC();
            if(NTmap.get(i).isGC()) GCcontent++;
        }
    }


    // Returns true if the Cigar code contains an 'N'. This indicates that the read contains a split.
    private boolean isSplitRead(Cigar cigar) {
        boolean ret = false;
        for(CigarElement ce : cigar.getCigarElements()) {
            if(ce.toString().contains("N")) {
                ret = true;
                break;
            }
        }
        return(ret);
    }


    public void addReadData(String rn, long readStart, long readEnd, byte[] readBases) {

        int j = 0;
        for(long i = readStart; i <= readEnd; i++) {
            if( NTmap.containsKey(i) ) {
                int nt = readBases[j]; // Bytes 65:A, 67:C, 71:G, 84:T
                NTmap.get(i).add(rn, nt);
            }
        }
    }


    public void addReadDataByStart(String rn, long readStart, int nt) {

        if(NTmap.containsKey(readStart)) {
            // Bytes 65:A, 67:C, 71:G, 84:T
            NTmap.get(readStart).add(rn, nt);
        }
    }


    public void printReadDepth() {
        for(Long i : NTmap.keySet()) {
            Site pos = NTmap.get(i);
            System.out.println(i + "\t" + pos.getReadDepth());
        }
    }




}
