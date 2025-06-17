import yagmail
from loguru import logger


def send_email_notification(cfg):
    """
    Send an email notification after the PDB extraction job completes.

    Args:
        new_structures (int): The number of new structures added.
        email_user (str): The sender's email address.
        email_password (str): The sender's email password.
        receiver_email (str): The receiver's email address.
        log_file (str): Path to the log file for logging email notifications.

    Returns:
        None
    """
    email_user = cfg.email.user
    email_password = cfg.email.password
    receiver_email = cfg.email.recipient

    subject = "PDB Extraction Job Completed"
    body = "The job has completed successfully."

    try:
        yag = yagmail.SMTP(email_user, email_password)
        yag.send(to=receiver_email, subject=subject, contents=body)
        logger.info(
            f"Email sent to {receiver_email} from {email_user} with subject '{subject}' and body '{body}'"
        )
    except Exception as e:
        logger.error(f"Failed to send email to {receiver_email}: {e}")
