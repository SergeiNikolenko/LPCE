import yagmail
from loguru import logger


def send_email_notification(
    new_structures: int,
    email_user: str,
    email_password: str,
    receiver_email: str,
    log_file: str,
):
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
    # Add a log file for the email notification process
    logger.add(log_file, format="{time} | {level} | {message}", level="INFO")

    subject = "PDB Extraction Job Completed"
    body = (
        f"The job has completed successfully. {new_structures} new ligands were added."
    )

    try:
        yag = yagmail.SMTP(email_user, email_password)
        yag.send(to=receiver_email, subject=subject, contents=body)
        logger.info(
            f"Email sent to {receiver_email} from {email_user} with subject '{subject}' and body '{body}'"
        )
    except Exception as e:
        logger.error(f"Failed to send email to {receiver_email}: {e}")
