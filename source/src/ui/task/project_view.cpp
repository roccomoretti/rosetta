// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/project_view.cpp
/// @brief  Project view class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#include "project_view.h"
#include "ui_project_view.h"

#include <ui/task/task.h>
#include <ui/task/task_view.h>
#include <ui/task/task_submit.h>
#include <ui/task/task_cancel_dialog.h>

#include <ui/config/config_dialog.h>

#include <QFileDialog>
#include <QDebug>

#include <set>

namespace ui {
namespace task {


ProjectView::ProjectView(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ProjectView)
{
    ui->setupUi(this);

	project_ = std::unique_ptr<Project>(new Project);  //std::make_shared<Project>();

	ui->tasks->setModel( project_->model() );
	//ui->tasks->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	//ui->tasks->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
	//ui->tasks->horizontalHeader()->setStretchLastSection(false);
	//ui->tasks->horizontalHeader()->resizeSection(1, 100);

	ui->tasks->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
	ui->tasks->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Fixed);

	// QItemSelectionModel *selectionModel = new QItemSelectionModel(model);
	// ui->project->setSelectionModel(selectionModel);
	// ui->project->selectionModel()->setCurrentIndex(ui->project->currentIndex(), QItemSelectionModel::SelectCurrent);
}

ProjectView::~ProjectView()
{
    delete ui;
}


void ProjectView::on_new_task_clicked()
{
    qDebug() << "ProjectView::on_new_task_clicked()";

	if( check_submit_requirements(*project_) ) {

		//QString label = QString("Task-%1").arg( project_->size() );

		QString label = QString("Task-%1").arg( project_->tasks().size() );
		auto task = std::make_shared<Task>();
		task->name(label);

		project_->add_task(task);

		auto ts = new TaskSubmit(task);
		ts->show();
	}

	// if( check_submit_requirements(*project_) ) {
	// 	static int task_counter = 0;
	// 	//QString label = QString("Task-%1").arg( project_->size() );
	// 	QString label = QString("Task-%1").arg( task_counter++ );
	// 	project_->add(label, std::make_shared<Task>(label) );
	// 	project_model_->set(project_);
	// }
}


void ProjectView::on_clone_task_clicked()
{
    qDebug() << "ProjectView::on_clone_task_clicked()";

	QItemSelectionModel *select = ui->tasks->selectionModel();

	if( select->hasSelection() ) {
		auto rows = select->selectedIndexes();

		std::set<int> selection;
		for(auto const & index : rows) selection.insert( index.row() );

		std::vector<TaskSP> tasks_to_clone;
		for(auto r : selection) tasks_to_clone.push_back( project_->task(r) );

		for(auto const & t : tasks_to_clone) {
			auto task = t->clone();

			project_->add_task(task);

			auto ts = new TaskSubmit(task);
			ts->show();
		}
	}
}

std::vector<TaskSP> ProjectView::get_selection()
{
	std::vector<TaskSP> r;

	QItemSelectionModel *select = ui->tasks->selectionModel();
	if( select->hasSelection() ) {
		auto rows = select->selectedIndexes();

		std::set<int> selection;
		for(auto const & index : rows) selection.insert( index.row() );

		for(auto row : selection) r.push_back( project_->task(row) );
	}

	return r;
}


void ProjectView::on_cancel_task_clicked()
{
    qDebug() << "ProjectView::on_cancel_task_clicked()";

	for(auto const & task : get_selection()) {
		TaskCancelDialog tcd;
		tcd.run(task);
	}
}


void ProjectView::on_delete_task_clicked()
{
    qDebug() << "ProjectView::on_delete_task_clicked()";

	auto tasks_to_delete = get_selection();
	for(auto const & t : tasks_to_delete) project_->erase(t);


	// QItemSelectionModel *select = ui->tasks->selectionModel();

	// if( select->hasSelection() ) {
	// 	auto rows = select->selectedIndexes();

	// 	std::set<int> selection;
	// 	for(auto const & index : rows) selection.insert( index.row() );

	// 	std::vector<TaskSP> tasks_to_delete;
	// 	for(auto r : selection) {
	// 		qDebug() << "row:" << r;
	// 		tasks_to_delete.push_back( project_->task(r) );
	// 	}
	// 	for(auto const & t : tasks_to_delete) project_->erase(t);
	// }
}


void ProjectView::on_action_new_project_triggered()
{
	project_ = std::unique_ptr<Project>(new Project);
	ui->tasks->setModel( project_->model() );
	raise();
}


void ProjectView::on_action_open_project_triggered()
{
    qDebug() << "ProjectView::on_project_open()";

	QString file_name = QFileDialog::getOpenFileName(this, tr("Open Project"), "", tr("RosettaUI Projects (*.rosetta)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_name.isEmpty() ) {
		QFile file(file_name);

		if (!file.open(QIODevice::ReadOnly) ) return;

		QDataStream in( &file );

		project_ = std::unique_ptr<Project>(new Project);
		in >> *project_;
		file.close();

		ui->tasks->setModel( project_->model() );

		project_->file_name(file_name);
		raise();
	}
}

void ProjectView::on_action_save_project_triggered()
{
    qDebug() << "ProjectView::on_project_save()";
	save_project(*project_, /* always_ask_for_file_name =*/ false);
}

void ProjectView::on_action_save_project_as_triggered()
{
    qDebug() << "ProjectView::on_project_save_as()";
	save_project(*project_, /* always_ask_for_file_name =*/ true);
}


void ProjectView::on_action_preferences_triggered()
{
    qDebug() << "ProjectView::on_action_preferences_triggered()";

	config::ConfigDialog * preferences = new config::ConfigDialog(this);
    preferences->show();
}

void ProjectView::on_tasks_doubleClicked(const QModelIndex &index)
{
    qDebug() << "ProjectView::on_project_doubleClicked";

	//QModelIndex index = ui->project_view->currentIndex();

	if( index.isValid() ) {
		TaskSP task = project_->task( index.row() );

		if(task) {
			if( task->state() == Task::State::none ) {
				auto ts = new TaskSubmit(task);
				ts->show();
			} else {
				auto tv = new TaskView(task);
				tv->show();
			}
		}
	}
}


} // namespace task
} // namespace ui
