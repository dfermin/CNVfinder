package umich.CNVfinder;

public class Coord {

    public String chrom;
    public long start;
    public long end;
    public long size;

    public Coord(String arg, int flank) {
        if(arg.startsWith("chr")) {
            String tmp = arg.replace("chr", "");
            arg = tmp;
        }

        chrom = arg.substring(0, arg.indexOf(":") );
        start = Long.valueOf( arg.substring(arg.indexOf(":")+1, arg.indexOf("-")) );
        end = Long.valueOf( arg.substring(arg.indexOf("-")+1) );

        if(flank == 0) {
            long s = end - start + 1;
            int local_flank = (int) (s / 2) + 1;
            start -= local_flank;
            end += local_flank;
        }
        else {
            start -= flank;
            end += flank;
        }

        size = end - start + 1;
    }
}
