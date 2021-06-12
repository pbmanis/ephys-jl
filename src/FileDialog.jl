module FileDialog
"""
From the web: 
https://discourse.julialang.org/t/qml-filedialog/10125
very slightly modified and reformatted

"""

using QML
export uigetfile, uigetdir, uisavefile


    function filedialog(;title::String="File Dialog",
        foldermode::Bool=false,
        multiselect::Bool=false,
        filter::Array{String}=["All files (*)"],
        folder::String=pwd(),
        savemode::Bool=false)

        qml_data = QByteArray("""
        import QtQuick 2.2
        import QtQuick.Dialogs 1.0
        import QtQuick.Controls 1.0
        import org.julialang 1.0

        ApplicationWindow {
          title: "FileDialog"
          width: 640
          height: 480
          visible: false

          FileDialog {
              id: fileDialog
              title: Title
              selectMultiple: MultiSelect
              selectFolder: FolderSelect
              selectExisting: SelectExisting
              nameFilters: Filter
              folder: Folder
              onAccepted: {

                  Julia.getfilelist(fileDialog.fileUrls)
                  Qt.quit()
              }
              onRejected: {
                  console.log("Canceled")
                  Qt.quit()
              }
              Component.onCompleted: visible = true
          }
        }
        """)

        filelist=AbstractArray{}
        function getfilelist(uri_list)
          filelist=uri_list
        end
        @qmlfunction getfilelist

        qengine = init_qmlengine()
        qcomp = QQmlComponent(qengine)
        set_data(qcomp, qml_data, "")
        set_context_property(qmlcontext(),"Title",title)
        set_context_property(qmlcontext(),"FolderSelect",foldermode)
        set_context_property(qmlcontext(),"MultiSelect",multiselect)
        set_context_property(qmlcontext(),"Filter",filter)
        set_context_property(qmlcontext(),"Folder",folder)
        set_context_property(qmlcontext(),"SelectExisting",!savemode)
        create(qcomp, qmlcontext());

        exec()
        return filelist
    end

    function uigetfile(;title::String="Select a File",
        multiselect::Bool=false,
        filter::Array{String}=["All files (*)"],
        folder::String=pwd())
        file=filedialog(;title=title,
            foldermode=false,
            multiselect=multiselect,
            filter=filter,
            folder=folder,
            savemode=false)
        return file
    end

    function uigetdir(;title::String="Select a Folder",
        folder::String=pwd())
        file=filedialog(;title=title,
            foldermode=true,
            multiselect=false,
            folder=folder,
            savemode=false)
        return file
    end

    function uisavefile(;title::String="Save File As",
        filter::Array{String}=["All files (*)"],
        folder::String=pwd())
        file=filedialog(;title=title,foldermode=false,
             multiselect=false,
             filter=filter,
             folder=folder,
             savemode=true)
        return file
    end

end