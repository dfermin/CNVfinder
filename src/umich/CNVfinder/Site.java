package umich.CNVfinder;

import java.util.HashSet;

class Site {
    private HashSet<String> reads;
    private int readDepth;
    private int GC;
    private boolean isGC;

    public Site() {
        reads = new HashSet<>();
        readDepth = reads.size();
        GC = 0; // keeps track of the number of reads that call this site a G or a C
        isGC = false;
    }

    public void add(String rn, int ntCode) {
        reads.add(rn);
        readDepth = reads.size();
        if(ntCode == 67 || ntCode == 71) GC++;
    }

    // determine if this site should be called as G/C based upon the number of reads supporting that call
    public void callGC() {
        double frac = ((double) GC) / ((double) readDepth);
        if(frac > 0.5) isGC = true; // if more than half of the reads say this site is G/C we declare it to be so.
    }

    public boolean isGC() { return isGC; }

    public int getReadDepth() { return readDepth; }

    public HashSet<String> getAllRead() { return reads; }

    public void setReadDepth(int rd) { readDepth = rd; }
    public void set_isGC(boolean gc) { isGC = gc; }
    public boolean get_isGC() { return isGC; }

}

