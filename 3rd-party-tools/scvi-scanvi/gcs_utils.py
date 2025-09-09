from google.cloud import storage

def get_bucket_blob(bucket_name, remote_path):
    """Downloads file from Google Cloud Storage bucket

    :param bucket_name: GCS bucket name
    :param file_path: local file path
    """
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_name)
    return bucket.blob(remote_path)

def download_from_bucket(bucket_name, remote_path, local_path):
    """Downloads file from Google Cloud Storage bucket

    :param bucket_name: GCS bucket name
    :param file_path: local file path
    """
    print(f'Downloading file: {remote_path}')
    blob = get_bucket_blob(bucket_name, remote_path)
    blob.download_to_filename(local_path)
    print(f"{remote_path} downloaded to {local_path}.")

def delocalize_file(bucket_name, file_to_delocalize, bucket_destination):
    """Writes local file to Google bucket

    :param bucket_name: GCS bucket name
    :param file_to_delocalize: local file to delocalize
    :param bucket_destination: path in Google bucket to write file to
    """
    print(f'Uploading file: {file_to_delocalize}')
    blob = get_bucket_blob(bucket_name, bucket_destination)
    blob.upload_from_filename(file_to_delocalize)
    print(f"File {file_to_delocalize} uploaded to {bucket_destination}.")

def get_bucket_and_path(remote_file):
    """Get GCS bucket name and path from remote file

    :param remote_file: remote file path
    """
    bucket = remote_file[5:].split("/")[0]
    path = remote_file[5:].split("/")[1]
    return bucket, path

def is_remote_file(file_path):
    """Check if file exists in GCS bucket"""
    is_gs_url = file_path[:5] == "gs://"
    if is_gs_url:
        bucket_name, remote_path = get_bucket_and_path(file_path)
        blob = get_bucket_blob(bucket_name, remote_path)
        return blob.exists()
    else:
        return False

def pull_all_files(file_list):
    """Pull all files from GCS bucket

    :param file_list: list of remote file paths
    """
    for f in file_list:
        if is_remote_file(f):
            bucket_name, remote_path = get_bucket_and_path(f)
            filename = f.split('/')[-1]
            download_from_bucket(bucket_name, remote_path, filename)
