import os

import yagmail
from config.settings import EMAIL_PASSWORD, EMAIL_USER, RECEIVER_EMAIL


def send_email_notification(new_structures: int):
    subject = "PDB Extraction Job Completed"
    body = (
        f"The job has completed successfully. {new_structures} new ligands were added."
    )

    yag = yagmail.SMTP(EMAIL_USER, EMAIL_PASSWORD)
    yag.send(to=RECEIVER_EMAIL, subject=subject, contents=body)

    print(
        f"Email sent to {RECEIVER_EMAIL} from {EMAIL_USER} with subject '{subject}' and body '{body}'"
    )
