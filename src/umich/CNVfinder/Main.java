package umich.CNVfinder;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class Main {

    public static HashMap<String, TreeMap<Long, Integer>> safidMap = null;
    public static HashMap<Double, ArrayList<Double>> all_percentGCmap = null;
    public static HashMap<Double, Double> percentGCmap = null;
    public static ArrayList<Double> allMedianRC = null;
    public static Coord targetRegionCoord = null;
    public static int windowSize = 100; // default
    public static int mapq_cutoff = 30;
    public static double overall_median_RC;
    public static boolean normalize = true;
    public static boolean verbose = false;


    public static void main(String[] args) throws IOException {

        if(args.length < 3) {
            //System.err.println("\nUSAGE: java -jar CNVfinder.jar -w <window_size (default: 100) [-n] -v -q <MAPQ cut off (default: 30)> -f <flank_length> -d <folder_of_bam_files> -c <chr:start-end>\n");
            System.err.println("\nUSAGE: java -jar CNVfinder.jar -i -c [-w -v -q -f]");
            System.err.println("\t-i   either a single BAM file or a folder of BAM files to process (REQUIRED)");
            System.err.println("\t-c   target region in CHROM:START-END format (REQUIRED)\n");
            System.err.println("\t-w   window size, default: 100");
            System.err.println("\t-q   MAPQ score cut off, default: 30");
            System.err.println("\t-f   flank size, +/- this amount will surround target region, default: flank by +/- half size of region");
            System.err.println("\n");
            System.exit(1);
        }

        String bamInput = null;
        String targetRegionStr = null;
        int flank = 0;

        for(int i = 0; i < args.length; i++) {
            if(args[i].equalsIgnoreCase("-i")) bamInput = args[++i];
            if(args[i].equalsIgnoreCase("-c")) targetRegionStr = args[++i];
            if(args[i].equalsIgnoreCase("-w")) windowSize = Integer.valueOf( args[++i] );
            if(args[i].equalsIgnoreCase("-f")) flank = Integer.valueOf(args[++i]);
            if(args[i].equalsIgnoreCase("-q")) mapq_cutoff = Integer.valueOf(args[++i]);
            if(args[i].equalsIgnoreCase("-n")) normalize = false;
            if(args[i].equalsIgnoreCase("-v")) verbose = true;
        }

        System.err.println("\nRead MAPQ must be >= " + mapq_cutoff + "\n");

        targetRegionCoord = new Coord(targetRegionStr, flank);
        recordTargetWindows(bamInput, 1);

        if(normalize) {
            System.err.println("Adjusting for GC content...");
            calcGCnormalizationValues();
            recordTargetWindows(bamInput, 2);
        }

        printReadDepth();
    }


    // Function computes the median read counts for all windows having the same %GC and the median read count across all windows
    private static void calcGCnormalizationValues() {

        overall_median_RC = getMedian(allMedianRC);

        for(double k : all_percentGCmap.keySet()) {
            percentGCmap.put(k, getMedian(all_percentGCmap.get(k)));
        }
    }


    public static double getMedian(ArrayList<Double> a) {
        double ret = 0;

        Collections.sort(a); // you need 'a' to be sorted low-to-high

        if(a.size() > 0) {
            int middle = a.size() / 2;

            if(a.size() % 2 == 1) {
                ret = a.get(middle);
            }
            else {
                ret = (a.get((middle-1)) + a.get(middle)) / 2.0;
            }
        }
        return(ret);
    }


    private static void printReadDepth() {

        // Get SAFIDs out. You need to keep the order constant
        int N = safidMap.size();
        String[] safids = new String[N];
        int i = 0;
        for(String s : safidMap.keySet()) { safids[i++] = s; }

        // Print header
        System.out.println("coord\t" + String.join("\t", safids));

        long stop = (targetRegionCoord.end - windowSize);
        for(long pos = targetRegionCoord.start; pos <= stop; pos += windowSize) {
            for (i = 0; i < N; i++) {

                TreeMap<Long, Integer> RD = safidMap.get( safids[i] );
                if(i == 0) {
                    System.out.print(targetRegionCoord.chrom + ":" + pos + "\t");
                }

                System.out.print(RD.get(pos));
                if(i < N) System.out.print("\t");
            }
            System.out.print("\n");
        }
    }


    private static void recordTargetWindows(String bamInput, int iteration) throws IOException {

        File bamSrc = new File(bamInput);
        if( !bamSrc.exists() ) {
            System.err.println("\nERROR: Unable to access " + bamInput + "\n");
            System.exit(2);
        }

        safidMap = new HashMap<>();

        // Initialize all static variables here
        // You only have to do this on the first iteration.
        // If you do %GC normalization this data will not be needed on the second iteration
        if(iteration == 1) {
            all_percentGCmap = new HashMap<>();
            percentGCmap = new HashMap<>();
            allMedianRC = new ArrayList<>();
            double[] bins = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
            for (double bin : bins) {
                all_percentGCmap.put(bin, new ArrayList<>());
                percentGCmap.put(bin, 0.0);
            }
        }


        if( bamSrc.isDirectory() ) {
            // Iterate over each BAM file recording the read evidence for each NT position of the target window.
            // We store the target window for each sample into the 'safidMap'
            for (String curBamFileName : bamSrc.list()) {
                if (curBamFileName.endsWith(".bam")) {
                    if ((iteration == 1) && verbose) System.err.println(curBamFileName);

                    String k = curBamFileName.split("\\.")[0];

                    String fp = bamSrc + "/" + curBamFileName;

                    Window curWindow = new Window(targetRegionCoord);
                    curWindow.init();
                    curWindow.parseBAM(fp);
                    curWindow.getReadDepthByWindow(windowSize);


                    if (iteration == 1) {
                        safidMap.put(k, curWindow.getRDmap());
                        recordGC(curWindow.getPercentGC(), curWindow.getMedianReadCount());
                        allMedianRC.add(curWindow.getMedianReadCount());
                    } else {
                        curWindow.adjustReadCounts();
                        safidMap.put(k, curWindow.getRDmap());
                    }
                }
            }
            if (iteration == 1) System.err.println("Read in " + safidMap.size() + " bam files.");
        }
        else { // user provided a single BAM file

            System.err.println(bamSrc.getName());
            String k = bamSrc.getName().split("\\.")[0];
            Window curWindow = new Window(targetRegionCoord);
            curWindow.init();
            curWindow.parseBAM(bamInput);
            curWindow.getReadDepthByWindow(windowSize);

            // iteration 1
            safidMap.put(k, curWindow.getRDmap());
            recordGC(curWindow.getPercentGC(), curWindow.getMedianReadCount());
            allMedianRC.add(curWindow.getMedianReadCount());

            //iteration 2
            curWindow.adjustReadCounts();
            safidMap.put(k, curWindow.getRDmap());
        }
    }



    private static void recordGC(double percentGC, double medianRC) {

        double[] bins = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
        for(int j = 1; j < bins.length - 1; j++) {
            int i = j - 1;
            if( (percentGC >= bins[i]) && (percentGC < bins[j]) ) {
                all_percentGCmap.get(bins[i]).add(medianRC);
            }
        }
    }
}
