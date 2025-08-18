import os
import time
import pandas as pd
import pysam
from collections import defaultdict
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email import encoders

# CONFIG
WATCH_FOLDER = "/path/to/sequencing/runs"
CHECK_INTERVAL = 24 * 60 * 60  # 1 day in seconds
OUTPUT_FILE = "variant_report.csv"

# Email Settings
SMTP_SERVER = "smtp.gmail.com"
SMTP_PORT = 587
SENDER_EMAIL = "your_email@gmail.com"
SENDER_PASSWORD = "your_app_password"  # use app password for Gmail
RECIPIENT_EMAILS = ["recipient1@email.com", "recipient2@email.com"]

# Thresholds
MIN_READS = 5
MIN_COVERAGE = 20
VAF_TOLERANCE = 0.02  # 2%

def send_email(report_file):
    """Send the generated report as an email attachment."""
    msg = MIMEMultipart()
    msg['From'] = SENDER_EMAIL
    msg['To'] = ", ".join(RECIPIENT_EMAILS)
    msg['Subject'] = "Daily Variant Report"

    # Email body
    body = "Attached is the latest sequencing variant report."
    msg.attach(MIMEText(body, 'plain'))

    # Attach CSV
    with open(report_file, "rb") as f:
        part = MIMEBase('application', 'octet-stream')
        part.set_payload(f.read())
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', f'attachment; filename={os.path.basename(report_file)}')
        msg.attach(part)

    # Send
    try:
        server = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
        server.starttls()
        server.login(SENDER_EMAIL, SENDER_PASSWORD)
        server.sendmail(SENDER_EMAIL, RECIPIENT_EMAILS, msg.as_string())
        server.quit()
        print(f"Email sent successfully to {RECIPIENT_EMAILS}")
    except Exception as e:
        print(f"Error sending email: {e}")


def vcf_to_variants(vcf_file, sample_name):
    """Extract variants from VCF file with filters applied."""
    variants = []
    vcf = pysam.VariantFile(vcf_file)

    for record in vcf.fetch():
        for sample in record.samples:
            data = record.samples[sample]
            depth = data.get("DP", 0)
            alt_count = data.get("AD", [0, 0])[1] if "AD" in data else 0
            vaf = alt_count / depth if depth > 0 else 0.0
            qual = record.qual if record.qual else 0

            vartype = "SNP" if len(record.ref) == 1 and len(record.alts[0]) == 1 else "INDEL"

            if (alt_count >= MIN_READS and depth >= MIN_COVERAGE and qual >= 20):
                variants.append((f"{record.chrom}:{record.pos}{record.ref}>{record.alts[0]}",
                                 vartype, vaf, qual, depth, alt_count, sample_name))
    return variants


def process_vcfs(folder):
    """Process all VCFs in a folder and generate a report of shared variants within 2% VAF."""
    all_variants = defaultdict(list)

    for file in os.listdir(folder):
        if file.endswith(".vcf") or file.endswith(".vcf.gz"):
            sample_name = file.split(".")[0]
            vcf_path = os.path.join(folder, file)
            variants = vcf_to_variants(vcf_path, sample_name)
            for var in variants:
                variant_id, vartype, vaf, qual, depth, alt_count, sample = var
                all_variants[(variant_id, vartype)].append((vaf, sample))

    report_rows = []
    for (variant_id, vartype), values in all_variants.items():
        if len(values) > 1:
            vafs = [v[0] for v in values]
            if max(vafs) - min(vafs) <= VAF_TOLERANCE:
                samples = [v[1] for v in values]
                report_rows.append({
                    "Variant": variant_id,
                    "Type": vartype,
                    "VAF Range": f"{min(vafs):.3f}-{max(vafs):.3f}",
                    "Samples": ",".join(samples)
                })

    df = pd.DataFrame(report_rows)
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"Report saved: {OUTPUT_FILE}")

    # Send email notification with report
    if not df.empty:
        send_email(OUTPUT_FILE)
    else:
        print("No variants found to report. Email not sent.")


def monitor_folder():
    """Check folder once a day for new sequencing runs and process them."""
    seen_runs = set()

    while True:
        for run in os.listdir(WATCH_FOLDER):
            run_path = os.path.join(WATCH_FOLDER, run)
            if os.path.isdir(run_path) and run not in seen_runs:
                vcf_files = [f for f in os.listdir(run_path) if f.endswith(".vcf") or f.endswith(".vcf.gz")]
                if vcf_files:
                    print(f"Processing run: {run}")
                    process_vcfs(run_path)
                    seen_runs.add(run)
        time.sleep(CHECK_INTERVAL)


if __name__ == "__main__":
    monitor_folder()
