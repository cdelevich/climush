import ftplib
from abc import ABC, abstractmethod
from typing import override


class OpenData(ABC):

    @property
    @abstractmethod
    def host_server(self)->str:
        """URL of the server to which data will be sent.

        Returns:
            A string of the host server's URL.

        """
        pass

class TransferMethod(ABC):
    """Abstract class on which all data transfer methods are based."""

    @property
    @abstractmethod
    def requires_auth(self)->bool:
        """Does this transfer method require authentication?

        If transferring files from client (climush user) to host (data
        repository), client authentication will be required.
        """
        pass

class FileTransferProtocol(TransferMethod):
    """File Transfer Protocol (FTP) data transfers"""

    def __init__(
            self,
            host_url: str,
            username: str,
            password: str,
            account: str,
            timeout: float|None,
            src_address: tuple|None,
            encoding,
    ):
        """

        :param host_url:
        :param username:
        :param password:
        :param account:
        :param timeout:
        :param src_address:
        :param encoding:
        """

    @abstractmethod
    def ftp(self):
        init_ftp = FTP(
            host='',
            user='',
            passwd='',
            acct='',
            timeout=None,
            source_address=None,
            encoding='utf-8',
        )
        return init_ftp


class NCBI(OpenData, FileTransferProtocols):
    """Upload sequencing data to NCBI.

    NCBI uses File Transfer Protocol (FTP) to accept files from clients.
    """

    @override
    def host_server(self):
        return 'ftp-private.ncbi.nlm.nih.gov'

    @override
    def ftp(self):