from base64 import b64decode, b64encode
from datetime import datetime
from hashlib import sha256
from hmac import new as new_hmac
from random import SystemRandom
from struct import pack, unpack

from Crypto.Cipher import AES

ENCRYPTION_KEY = b"*EnGiNeKey_4224_yeKeNiGnE*.D0N0tGuess///e"


def get_hmac(input: bytes) -> str:
    """
    Return a SHA256-based hmac
    """
    return new_hmac(ENCRYPTION_KEY, input, sha256).hexdigest()


def unpad(paddedtext: bytes) -> bytes:
    """
    Unpad a BioNumerics-style padded text
    """
    if len(paddedtext) % AES.block_size:
        raise RuntimeError("Not BN encrypted text")
    if b"ENCR01" != paddedtext[:6]:
        raise RuntimeError("Not BN encrypted text")
    if len(paddedtext) < 10:
        raise RuntimeError("Not BN encrypted text")
    (cleartext_len,) = unpack("i", paddedtext[6:10])
    if len(paddedtext) < 10 + cleartext_len:
        raise RuntimeError("Not BN encrypted text")
    return paddedtext[10 : 10 + cleartext_len]


def decrypt(ciphertext: bytes, iv: bytes):
    """
    Decrypt a BioNumerics-style encrypted text
    """
    cipher = AES.new(ENCRYPTION_KEY[:32], AES.MODE_CBC, iv=iv)
    paddedtext = cipher.decrypt(ciphertext)
    return unpad(paddedtext)


def decode_and_decrypt(base64text: bytes) -> bytes:
    """
    Decode and decrypt a base64-encoded BioNumerics-style encrypted text
    """
    decoded = b64decode(base64text)
    iv = decoded[: AES.block_size]
    ciphertext = decoded[AES.block_size :]
    return decrypt(ciphertext, iv)


def pad(cleartext: bytes) -> bytes:
    """
    Pad a text with BioNuemrics-specific padding
    """
    cleartext_len = len(cleartext)
    mod = (10 + cleartext_len) % AES.block_size
    while mod:
        cleartext += cleartext[: AES.block_size - mod]
        mod = (10 + len(cleartext)) % AES.block_size
    return b"ENCR01" + pack("i", cleartext_len) + cleartext


def encrypt(cleartext: bytes, iv: bytes = None) -> bytes:
    """
    Encrypt a text with BioNumerics-specific encryption
    """
    cipher = AES.new(ENCRYPTION_KEY[:32], AES.MODE_CBC, iv=iv)
    paddedtext = pad(cleartext)
    return cipher.encrypt(paddedtext), cipher.iv


def encrypt_and_encode(cleartext: bytes) -> bytes:
    """
    BioNumerics-style encrypt and then Base-64 encode a text
    """
    ciphertext, iv = encrypt(cleartext)
    return b64encode(iv + ciphertext)


def sprinkle(input: str) -> str:
    """
    Sprinkel a string with random characters
    """
    positions = (3, 1, 3, 1, 5, 1, 2, 6, 5, 3, 2, 5)  # alsmost pi :)
    client_id = ""
    pos_index = 0
    prev_pos = 0
    while pos_index < len(positions) and prev_pos < len(input):
        next_pos = min(len(input), prev_pos + positions[pos_index])
        client_id += input[prev_pos:next_pos]  # take a piece of original str
        if next_pos == prev_pos + positions[pos_index]:
            client_id += str(
                SystemRandom().randrange(10)
            )  # add some random stuff at fixed positions
        prev_pos = next_pos
        pos_index += 1
    if prev_pos < len(input):
        client_id += input[prev_pos:]
    return client_id


def unsprinkle(input: str) -> str:
    """
    Unsprinkel random characters from a string
    """
    positions = (3, 1, 3, 1, 5, 1, 2, 6, 5, 3, 2, 5)  # alsmost pi :)
    client_id = ""
    pos_index = 0
    prev_pos = 0
    while pos_index < len(positions) and prev_pos < len(input):
        next_pos = min(len(input), prev_pos + positions[pos_index])
        client_id += input[prev_pos:next_pos]  # take a piece of original str
        if next_pos == prev_pos + positions[pos_index]:
            next_pos += 1
        prev_pos = next_pos
        pos_index += 1
    if prev_pos < len(input):
        client_id += input[prev_pos:]
    return client_id


def get_client_id(serial: str, password: str, timestamp: datetime = None) -> bytes:
    """
    Creates encoded login/pwd to send to the nomenclature service.
    """
    epoch = datetime(2010, 1, 1)
    now = timestamp or datetime.utcnow()
    seconds = (now - epoch).total_seconds()  # this is a float !
    input = serial + str(int(seconds)).zfill(11) + password
    client_id = sprinkle(input)
    return encrypt_and_encode(client_id.encode("utf-8"))
