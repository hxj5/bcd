# subset_xf_tag.py


import pysam
import sys



def subset_bam_by_xf_tag(in_sam_fn, out_sam_fn, xf_tag = "xf"):
    in_sam = pysam.AlignmentFile(in_sam_fn, "r")
    out_sam = pysam.AlignmentFile(out_sam_fn, "wb", template = in_sam)
    
    n_in, n_out = 0, 0
    for read in in_sam.fetch():
        n_in += 1
        if not read.has_tag(xf_tag):
            continue
        flag = read.get_tag(xf_tag)
        flag = int(flag)
        if flag in (17, 25):
            out_sam.write(read)
            n_out += 1
            
    in_sam.close()
    out_sam.close()
    return((n_in, n_out))
    


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: %s <in_bam> <out_bam>" % sys.argv[0])
        sys.exit(1)
        
    in_sam_fn, out_sam_fn = sys.argv[1:3]
    n_in, n_out = subset_bam_by_xf_tag(in_sam_fn, out_sam_fn)
    print("n_reads_in = %d; n_reads_out = %d." % (n_in, n_out))
