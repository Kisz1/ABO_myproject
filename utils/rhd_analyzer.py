from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
import streamlit as st

class RHDAnalyzer:
    def __init__(self, gb_path="utils/data/rhd_referance.gb", reference_seq=None):
        """
        Load reference sequence from a GenBank file or from a provided sequence.
        """
        self.reference_seq = None
        self.exon_map = []

        if reference_seq is not None:
            self.reference_seq = str(reference_seq)
            self.exon_map = []
        else:
            try:
                record = SeqIO.read(gb_path, "genbank")
                self.reference_seq = str(record.seq)
                # ดึงตำแหน่ง Exon อัตโนมัติจากไฟล์ GenBank
                self.exon_map = self._get_exon_map(record)
            except Exception as e:
                st.error(f"ไม่สามารถโหลดไฟล์ Reference ได้: {e}")

    def _get_exon_map(self, record):
        exons = []
        for feature in record.features:
            if feature.type == "exon":
                # เก็บเป็น tuple (start, end)
                start = int(feature.location.start)
                end = int(feature.location.end)
                exons.append((start, end))
        exons.sort() # เรียงจาก Exon 1 ไป 10
        return exons

    def map_genomic_to_coding(self, genomic_idx):
        """
        แปลงตำแหน่งจากสาย Genomic (0-58000) เป็น Coding Sequence (c.)
        """
        if not self.exon_map:
            return genomic_idx + 1

        coding_pos = 0
        for start, end in self.exon_map:
            if start <= genomic_idx < end:
                return coding_pos + (genomic_idx - start + 1)
            coding_pos += (end - start)
        return None

    def analyze(self, query_seq):
        """
        ฟังก์ชันหลัก: รับสายเบสของคนไข้แล้วคืนค่ารายการความผิดปกติและความเหมือน
        โดยตรวจทั้งสายปกติและ reverse complement แล้วเลือกผลลัพธ์ที่มี identity สูงสุด
        """
        forward_result = self._analyze_sequence(str(query_seq))
        reverse_result = self._analyze_sequence(str(Seq(query_seq).reverse_complement()))

        if reverse_result["identity"] > forward_result["identity"]:
            reverse_result["strand"] = "reverse"
            return reverse_result

        forward_result["strand"] = "forward"
        return forward_result

    def _analyze_sequence(self, query_seq):
        alignments = pairwise2.align.localms(self.reference_seq, query_seq, 2, -1, -0.5, -0.1)

        if not alignments:
            return {"variants": [], "identity": 0.0, "score": 0.0}

        best_match = alignments[0]
        variants = self.extract_variants(best_match)
        identity = self._compute_alignment_identity(best_match)
        score = best_match.score
        return {"variants": variants, "identity": identity, "score": score}

    def _compute_alignment_identity(self, alignment):
        """
        Compute identity % only for aligned bases (ignoring gaps).
        """
        ref_aligned = alignment[0]
        query_aligned = alignment[1]
        matches = 0
        aligned_positions = 0

        for ref_b, query_b in zip(ref_aligned, query_aligned):
            if ref_b == "-" or query_b == "-":
                continue
            aligned_positions += 1
            if ref_b == query_b:
                matches += 1

        if aligned_positions == 0:
            return 0.0

        return (matches / aligned_positions) * 100.0

    def extract_variants(self, alignment):
        """
        Extract variants only from positions where both ref and query are present (no gaps).
        """
        ref_aligned = alignment[0]
        query_aligned = alignment[1]
        start_in_ref = alignment[3]
        
        variants = []
        current_ref_idx = start_in_ref
        
        for i in range(len(ref_aligned)):
            ref_b = ref_aligned[i]
            que_b = query_aligned[i]
            
            if ref_b != "-":
                current_ref_idx_pos = current_ref_idx
                current_ref_idx += 1
            else:
                continue
            
            if que_b == "-":
                continue
            
            if ref_b != que_b:
                c_pos = self.map_genomic_to_coding(current_ref_idx_pos)
                if c_pos is None:
                    continue

                prefix = "c." if self.exon_map else "g."
                if que_b == "-":
                    variants.append(f"{prefix}{c_pos}del{ref_b}")
                else:
                    variants.append(f"{prefix}{c_pos}{ref_b}>{que_b}")
            
        return variants
    
# --- ส่วนสำหรับทดสอบรัน (วางไว้ท้ายไฟล์) ---
if __name__ == "__main__":
    import os

    # 1. ระบุพาธไฟล์ (แก้ให้ตรงกับชื่อไฟล์จริงของคุณ)
    test_fasta = "sample_for_global_alignment/TSN20251113-010-00152_ABO-Chon_20251115-BAN10_C02_C02_i37.1.fasta"
    gb_ref = "utils/data/rhd_referance.gb"

    if not os.path.exists(test_fasta):
        print(f"❌ ไม่พบไฟล์ทดสอบ: {test_fasta}")
    elif not os.path.exists(gb_ref):
        print(f"❌ ไม่พบไฟล์อ้างอิง: {gb_ref}")
    else:
        print("🚀 เริ่มการทดสอบวิเคราะห์ RHD...")
        # 2. สร้างตัววิเคราะห์
        analyzer = RHDAnalyzer(gb_path=gb_ref)
        
        # 3. รันการวิเคราะห์
        with open(test_fasta, "r") as f:
            found_variants = analyzer.analyze(f)
        
        # 4. แสดงผล
        print("\n=== ผลการวิเคราะห์ ===")
        if found_variants:
            print(f"พบจุดกลายพันธุ์ {len(found_variants)} ตำแหน่ง:")
            for v in found_variants:
                print(f" - {v}")
        else:
            print("✅ ไม่พบความผิดปกติในส่วน Exon (ลำดับเบสตรงกับอ้างอิง)")