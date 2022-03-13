function getCookie(name) {
    let cookieValue = null;
    if (document.cookie && document.cookie !== '') {
        const cookies = document.cookie.split(';');
        for (let i = 0; i < cookies.length; i++) {
            const cookie = cookies[i].trim();
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}

function data() {
    const csrftoken = getCookie('csrftoken');
    var form = new FormData()
    const value1 = document.getElementById('random_paulis').value.trim();
    const value2 = document.getElementById('number_qubits').value.trim();
    console.log("yes or no : ", document.getElementById('error_plot').checked)
    if (value1 != '' && value2 != '') {

        console.log("call")
        form.append('random_paulis', document.getElementById('random_paulis').value.trim())
        form.append('physical_qubits', document.getElementById('number_qubits').value.trim())
        form.append('error_plot', document.getElementById('error_plot').checked)
        fetch('http://127.0.0.1:8000/data', {
            method: 'POST',
            body: form,
            headers: {
                'X-CSRFToken': csrftoken
            },
            mode: 'same-origin',
        }).then(data => data.json())
        .then(data => {
            var result = data['result']
            document.getElementById('result_random_paulis').value = result['number_of_cliques']
            document.getElementById('matrix').innerHTML = result['text']
            document.getElementById('image_one').src = `http://127.0.0.1:8000/static${result['image_one']}`
            document.getElementById('image_two').src = `http://127.0.0.1:8000/static${result['image_two']}`
            document.getElementById('text_one').innerHTML = result['text_one']
            document.getElementById('text_two').innerHTML = result['text_two']
            document.getElementById('result_div').classList.remove('d-none')
            var data = ''
            result?.table?.map(val => {
                data += `
                        <tr>
                            <td>${val.val1}</td>
                            <td>${val.val2}</td>
                        </tr>
                    `
            })
            document.getElementById('table_body').innerHTML = data
            document.getElementById('error_div').innerHTML = '';
            console.log(data)
        })
    } else {
        document.getElementById('error_div').innerHTML = 'Please fill all the fields';
    }
}