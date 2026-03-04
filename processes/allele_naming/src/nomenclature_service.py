from dataclasses import dataclass
from xml.etree.ElementTree import Element, fromstring, tostring

from requests import PreparedRequest, Request, Session

from src.encryption import get_client_id, get_hmac


class NomenclatureServiceError(Exception):
    ...


@dataclass
class NomenclatureService:
    url: str = ""
    project: str = ""
    password: str = ""
    serial: str = ""  # BioNumerics serial number linked to the project

    def submit_alleles(self, alleles: Element) -> Element:
        return self._send_request("AddAlleles", alleles)

    def get_organisms(self) -> list[str]:
        organisms = self._send_request("GetOrganisms")
        return [organism.findtext("ID") for organism in organisms.findall("organism")]

    def _send_request(self, page: str, xml: Element = None) -> Element:
        """
        Send a request to the nomenclature service, optionally passing an XML to post

        :param str page: name of the ashx-page
        :param Element xml: xml to post (default None)
        :return Element node with the response of the request
        """
        if xml is None:
            request = self._prepare_get_request(page)
        else:
            request = self._prepare_post_request(page, xml)

        with Session() as session:
            response = session.send(request)
        response.raise_for_status()
        xml_response = fromstring(response.content.decode("utf-8"))
        if xml_response.tag == "Error":
            raise NomenclatureServiceError(
                "Nomenclature request was not successful: "
                + (
                    xml_response.findtext("Message")
                    or xml_response.text
                    or "Unknown error occurred."
                )
            )
        return xml_response

    def _get_page_url(self, page: str) -> str:
        """
        Get the url for a .ashx page
        """
        return self.url + ("" if self.url.endswith("/") else "/") + page + ".ashx"

    def _prepare_get_request(self, page: str) -> PreparedRequest:
        """
        Create the request to send to the nomenclature service
        """
        url = self._get_page_url(page)
        params = self._create_url_query()
        headers = {"User-Agent": "Nextflow"}
        request = Request("GET", url, params=params, headers=headers)
        return request.prepare()

    def _prepare_post_request(self, page: str, xml: Element) -> PreparedRequest:
        """
        Create the request to send to the nomenclature service
        """
        url = self._get_page_url(page)
        data = tostring(xml)
        params = self._create_url_query(data)
        headers = {"User-Agent": "Nextflow", "Content-Type": "text/xml"}
        request = Request("POST", url, params=params, data=data, headers=headers)
        return request.prepare()

    def _create_url_query(self, data: bytes = b"") -> dict[str, str]:
        """
        Create the URL query to send to the nomenclature service
        """
        client_id = get_client_id(self.serial, self.password)
        to_sign = client_id + b"||" + data
        return {
            "version": "2",
            "id": client_id.decode("utf-8"),
            "db": self.project,
            "hmac": get_hmac(to_sign),
        }
